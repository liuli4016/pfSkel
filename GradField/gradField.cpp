// ----
// ----  Computes the potential field for a volume
// ----  Input: volume file, dimensions: X, Y, Z, output file name
// ----	 Output: normalized potential field:
//		 1 vector for each point in the volume
//
// Last change: Thu May 15 15:20:38 EDT 2003 by Nicu D. Cornea
//
//

// #define TRACE

#include "gradField.h"
#include "critPts.h"

inline bool SmoothSurfaceVectors(
  unsigned char *f, int L, int M, int N, 
  ForceVector *force, int smSteps, int ng[26], 
  double *totalChange, double *minChange, double *maxChange, int *nrChanged);

inline bool CheckNormals(unsigned char *fc, int L, int M, int N, 
			 ForceVector *force, int ng[26]);

inline bool CriticalPointsAreClose(
				   unsigned char *flags, 
				   ForceVector *ForceField,
				   unsigned char *cpMap,
				   int *ng,
				   int L, int M, int N);

bool SetNeighborsArray(int ng[26], int L, int M, int N);

bool ComputeSurfaceNormals(unsigned char *fc, int L, int M, int N, 
			   ForceVector *force);


bool ComputeSubSurfaceVectors(ForceVector *force, unsigned char *fc, 
			      int L, int M, int N, int ng[26]);

bool NormalsAdaptiveAlg(ForceVector *force, 
			unsigned char *fc, unsigned char *flagsCopy, 
			int L, int M, int N, int ng[26], int mSmSteps, 
			int *nrSmStepsPerformed);

bool CritPtsAdaptiveAlg(ForceVector *force, 
			unsigned char *fc, unsigned char *flagsCopy, 
			int L, int M, int N, int ng[26], int mSmSteps, 
			bool useProvidedForces,
			int *nrSmStepsPerformed);

///////////////////////////////////////////////////////////////////////////////
// Function CalculateGradientField
//   Computes the gradient field generated by propagating the boundary gradient
// towards the center of the object.
///////////////////////////////////////////////////////////////////////////////
bool CalculateGradientField(
	int L, int M, int N, 	      // [in] size of volume
	unsigned char* f, 	      // [in] volume flags
	ForceVector *force,		      // [out] force field	
	int smSteps /*= 10*/,             // [in] number of smoothing steps
	Adaptive_Smoothing_Algorithm adaptiveAlgorithm /*= ASA_NONE*/
	                              // [in] adaptive algorithm to use:
) {

  long slsz, sz;
  unsigned char *fc;
  bool retVal = true;
  
  //
  // check volume padding
  //
  if(!CheckVolumePadding(f, L, M, N, 1)) {
    printf("** Error - Object touches bounding box. Abort.\n");
    return false;
  }


  slsz = L*M;		// slice size
  sz = slsz*N;

  
  int ng[26];
  // pre-compute offsets of all the 26 neighbors of a voxel in this volume
  SetNeighborsArray(ng, L, M, N);
  
  
  //
  // make a copy of the flags - because they will be changed during this 
  // procedure and we don't want to affect the calling program
  //
  if((fc = new unsigned char[L*M*N]) == NULL) {
    printf("Error allocating memory !!\n");
    exit(1);
    return false;
  }
  memcpy(fc, f, L*M*N*sizeof(unsigned char));

  // flag interior, boundary and exterior voxels
  FlagVolume(fc, L, M, N);

  // set all forces to 0
  memset(force, 0, sizeof(ForceVector)*L*M*N);

  retVal = true;

  switch(adaptiveAlgorithm) {
  case ASA_NONE:
    {
      // compute surface normals
      ComputeSurfaceNormals(fc, L, M, N, force);
      
      // smooth using smSteps steps
      printf("Performing mandatory smoothing steps ...\n");

      double totalChange, minChange, maxChange;// records the amount of change 
      int nrChanged;
      totalChange = 0.00;
      minChange = 0.00;
      maxChange = 0.00;
      nrChanged = 0;
      if(!SmoothSurfaceVectors(fc, L, M, N, force, smSteps, ng, 
			       &totalChange, &minChange, &maxChange, 
			       &nrChanged))
      {
	retVal = false;
	break;
      }
      
#ifdef _DEBUG
      PrintElapsedTime("\tPF-3: computing the gradient field for surface voxels.");
#endif
      
      // compute sub-surface vectors
      ComputeSubSurfaceVectors(force, fc, L, M, N, ng);
      
      printf("Done.\n");
    }
    break;
  case ASA_NORMALS:
    {
      int nrSmStepsPerformed = 0;
      if(!NormalsAdaptiveAlg(force, fc, f, L, M, N, ng, smSteps, 
			     &nrSmStepsPerformed)) 
      {
	printf("NormalsAdaptiveAlg(...) failed !\n");
	retVal = false;
      }
    }
      break;
  case ASA_CRITICAL_POINTS:
    {
      int nrSmStepsPerformed = 0;
      if(!CritPtsAdaptiveAlg(force, fc, f, L, M, N, ng, smSteps,
			     false,
			     &nrSmStepsPerformed)) 
      {
	printf("CritPtsAdaptiveAlg(...) failed !\n");
	retVal = false;
      }
    }
    break;
  case ASA_NORMALS_AND_CRITICAL_POINTS:
    {
      int nrSmStepsPerformed = 0;
      int tmp;
      if(!NormalsAdaptiveAlg(force, fc, f, L, M, N, ng, smSteps, 
			     &nrSmStepsPerformed)) 
      {
	printf("NormalsAdaptiveAlg(...) failed !\n");
	retVal = false;
      }
      // after this is done, do the critical point algorithm
      // and actually use the provided forces before re-smoothing
      if(!CritPtsAdaptiveAlg(force, fc, f, L, M, N, ng, nrSmStepsPerformed, 
			     true,
			     &tmp)) 
      {
	printf("CritPtsAdaptiveAlg(...) failed !\n");
	retVal = false;
      }
    }
    break;
  } // switch

  printf("\n");
  
  delete [] fc;
    
#ifdef _DEBUG
  PrintElapsedTime("\tPF-3: computing the entire gradient field.");
#endif
  
  return retVal;
}


bool NormalsAdaptiveAlg(ForceVector *force, 
			unsigned char *fc, unsigned char *flagsCopy, 
			int L, int M, int N, int ng[26], int mSmSteps, 
			int *nrSmStepsPerformed)
{
  int nrSurfaceVoxels;
  int i, idx, sz;
  double totalChange, minChange, maxChange;// records the amount of change 
  int nrChanged;
  double prevMaxChanges[3];
  bool done;
  int steps;
  
  sz = L*M*N;
  *nrSmStepsPerformed = 0;

  // compute nr of surface voxels
  nrSurfaceVoxels = 0;
  for(idx=0; idx < sz; idx++) {
    if(fc[idx] == SURF) nrSurfaceVoxels++;
  }

  // set all forces to 0
  memset(force, 0, sizeof(ForceVector)*L*M*N);

  // compute surface normals
  ComputeSurfaceNormals(fc, L, M, N, force);
  
  //
  // do at least <mSmSteps> steps of smoothing
  // to cancel out noise from normal estimation
  printf("Performing mandatory smoothing steps ...\n");
  totalChange = 0.00;
  minChange = 0.00;
  maxChange = 0.00;
  nrChanged = 0;
  if(!SmoothSurfaceVectors(fc, L, M, N, force, mSmSteps, ng, 
			   &totalChange, &minChange, &maxChange, &nrChanged)) {
    return false;
  }
  printf("Done.\n");

  (*nrSmStepsPerformed) = mSmSteps;
  
  // check the normals
  if(CheckNormals(fc, L, M, N, force, ng)) {
    // we are done
    printf("No 2 neighboring normals are the same. ");
  }
  else {
  
    //
    // perform adaptive smoothing steps
    //
    prevMaxChanges[0] = 0.00; prevMaxChanges[1] = 0.00; 
    prevMaxChanges[2] = 0.00;
    
    i = 0;
    done = false;
    steps = 1;
    while(!done) {
      printf("Step: %d. Performing another %d steps of smoothing... \n", 
	     i, steps);
      i++;

      // smooth with another <steps> steps
      totalChange = 0.00;
      minChange = 0.00;
      maxChange = 0.00;
      nrChanged = 0;
      SmoothSurfaceVectors(fc, L, M, N, force, steps, ng, 
			   &totalChange, &minChange, &maxChange, &nrChanged);
      printf("(change: (%.3lf - %.7lf), tot: %.6lf, #: %d)\n", 
	     minChange, maxChange, totalChange, nrChanged);
      fflush(stdout);
      
      (*nrSmStepsPerformed) = (*nrSmStepsPerformed)  + steps;
      
      if(!CheckNormals(fc, L, M, N, force, ng)) {
	
	// terminating condition
	// ---
	// this loop can go on forever if vectors change back and forth 
	// from one step to the next
	// setting a threshold on the total change only may amount to too many
	// steps and be very slow.
	// 
	// We terminate when all the surface normals are changing at each step 
	// and the maximum change in a normal direction is below
	// 1 - abs(cos(0.1)) = 1.5230867123011970987520742872343e-6
	// that is, the maximum change is less than 0.1 degrees.
	// ----
	
	if((nrChanged == nrSurfaceVoxels) && (maxChange < 1.5231e-6)) {
	  printf("Terminating condition fulfilled.\n");
	  done = true;
	  break;
	}

	// double the number of steps 
	// steps = steps * 2;
      }
      else {
	done= true;
	printf("No 2 neighboring normals are the same.\n");
	break;
      }
      // printf("done.\n");
      
      //
      // stop after 500 steps
      if((*nrSmStepsPerformed) > 500) {
	// stop here - too many steps
	printf("Stopping after 500 steps.\n");
	done = true;
	break;
      }
    }
  }

#ifdef _DEBUG
  PrintElapsedTime("\tPF-3: computing the gradient field for surface voxels.");
#endif

  // compute sub-surface vectors
  ComputeSubSurfaceVectors(force, fc, L, M, N, ng);
  
  printf("\n");
 
  return true;
}




bool CritPtsAdaptiveAlg(ForceVector *force, 
			unsigned char *fc, unsigned char *flagsCopy, 
			int L, int M, int N, int ng[26], int mSmSteps, 
			bool useProvidedForces,
			int *nrSmStepsPerformed)
{
  double totalChange, minChange, maxChange;// records the amount of change 
  int nrChanged;
  int i,j, dx, dy, dz;
  DynamicArray<CriticalPoint> CritPts;
  
  bool done = false;
  bool zeroEigenvalue = false;
  bool cpAreClose = false;
  bool retVal = true;
  
  ForceVector *prevForces;

  retVal = true;
  *nrSmStepsPerformed = mSmSteps;
  
  if((prevForces = new ForceVector[L*M*N]) == NULL) {
    printf("Out of memory !\n");
    return false;
  }
  memset(prevForces, 0, sizeof(ForceVector)*L*M*N);
  
  
  // if useProvidedForces is true, then we can use the forces given as 
  // parameter and continue smoothing from there
  // if not, we need to compute the forces first
  if(!useProvidedForces) {
    // compute forces first
    // set all forces to 0
    memset(force, 0, sizeof(ForceVector)*L*M*N);
    
    // Compute force vectors with current smoothing
    // assume flags are ok
    ComputeSurfaceNormals(fc, L, M, N, force);
    // smooth with <steps> steps
    totalChange = 0.00;
    minChange = 0.00;
    maxChange = 0.00;
    nrChanged = 0;
    if(!SmoothSurfaceVectors(fc, L, M, N, force, mSmSteps, ng, 
			     &totalChange, &minChange, &maxChange, &nrChanged))
    {
      retVal = false;
      done = true;
    }
    *nrSmStepsPerformed = mSmSteps;
    // compute sub-surface vectors
    ComputeSubSurfaceVectors(force, fc, L, M, N, ng);
  }
  
  // copy forces in prevForces
  memcpy(prevForces, force, sizeof(ForceVector)*L*M*N);
  
  
  if(retVal != false) {
    while(!done) {
      
      //
      // compute critical points
      //
      CritPts.Reset();
      if(!GetCriticalPoints(force, L, M, N, flagsCopy, CritPts, 
			    /*inOut =*/ false, /*notCloseToSurf =*/ false,
			    /*maxRecursionDepth =*/ 10,
			    /*useNewton =*/ true))
      {
	printf("GetCriticalPoints(...) failed !\n");
	done = true;
	retVal = false;
	break;
      }
      /*
      // print critical points
      for(i=0; i < CritPts.GetNrElem(); i++) {
      printf("%.16lf %.16lf %.16lf %d\n",
      CritPts[i].position.x, CritPts[i].position.y, 
      CritPts[i].position.z, CritPts[i].type);
      }
      */
      // check the critical points
      
      // if any of the critical points has a 0 eigenvalue - not done
      zeroEigenvalue = false;
      for(i=0; (i < CritPts.GetNrElem()) && (!zeroEigenvalue); i++) {
	for(j=0; j < 3; j++) {
	  if(CritPts[i].eval[j] == 0) {
	    zeroEigenvalue = true;
	    break;
	  }
	}
      }
      
      if(!zeroEigenvalue) {
	// check that every pair of attracting nodes is at least 2 voxels away
	cpAreClose = false;
	for(i=0; (i < CritPts.GetNrElem()) && (!cpAreClose); i++) {
	  if(CritPts[i].type == CPT_ATTRACTING_NODE) {
	    for(j=i+1; (j < CritPts.GetNrElem()) && (!cpAreClose); j++) {
	      if(CritPts[j].type == CPT_ATTRACTING_NODE) {
		// check distance
		dx = 
		  abs(int(CritPts[i].position.x) - int(CritPts[j].position.x));
		dy = 
		  abs(int(CritPts[i].position.y) - int(CritPts[j].position.y));
		dz = 
		  abs(int(CritPts[i].position.z) - int(CritPts[j].position.z));
		
		if((dx < 2) && (dy < 2) && (dz < 2)) {
		  cpAreClose = true;
		  break;
		}
	      }
	    }
	  }
	}
      }
      
      if(zeroEigenvalue || cpAreClose || (CritPts.GetNrElem() == 0)) {
	// keep going
	//
	// stop after 500 steps
	if((*nrSmStepsPerformed) > 500) {
	  // stop here - too many steps
	  printf("Stopping after 500 steps.\n");
	  done = true;
	  break;
	}
      }
      else {
	// we are done
	done = true;
	break;
      } 
      
      // if we are not done, smooth again
      if(!done) {
	// smooth
	printf("Using %d smoothing steps ...\n", (*nrSmStepsPerformed)+1);
	
	// re-copy flags 
	memcpy(fc, flagsCopy, L*M*N*sizeof(unsigned char));
	
	// flag interior, boundary and exterior voxels
	FlagVolume(fc, L, M, N);
	
	// copy previous forces
	memcpy(force, prevForces, sizeof(ForceVector)*L*M*N);
	
	// smooth with 1 more step
	totalChange = 0.00;
	minChange = 0.00;
	maxChange = 0.00;
	nrChanged = 0;
	if(!SmoothSurfaceVectors(fc, L, M, N, force, 1, ng, 
				 &totalChange, &minChange, &maxChange, 
				 &nrChanged))
	{
	  retVal = false;
	  done = true;
	  break;
	}
	*nrSmStepsPerformed = (*nrSmStepsPerformed) + 1;
	// compute sub-surface vectors
	ComputeSubSurfaceVectors(force, fc, L, M, N, ng);
	
	// save new forces
	memcpy(prevForces, force, sizeof(ForceVector)*L*M*N);
      }
    }
  }
  
  delete [] prevForces;
  
  printf("Done after %d smoothing steps.\n", (*nrSmStepsPerformed));
  return retVal;
}




inline bool SmoothSurfaceVectors(
  unsigned char *f, int L, int M, int N, 
  ForceVector *force, int smSteps, int ng[26], 
  double *totalChange, double *minChange, double *maxChange, int *nrChanged)
{
  int Lm1, Mm1, Nm1;
  int i,j,k, bk, s, p;
  long idx, iidx, bidx, slsz, sz;
  double oxd, oyd, ozd;
  bool changed;
  double change = 0;

  slsz = L*M;
  sz = slsz*N;
  
  Lm1 = L - 1;
  Mm1 = M - 1;
  Nm1 = N - 1;

  *totalChange = 0;
  *minChange = 1;
  *maxChange = 0;
  *nrChanged = 0;

  //
  // first step: normalize all surface vectors
  //
  for(idx=0; idx < sz; idx++) {
    if(f[idx] == SURF) {
      Normalize(force + idx);
    }
  }

  
  // keep a copy of the original data (3 slices) needed for smoothing
  ForceVector *buffer = NULL;
  if((buffer = new ForceVector[slsz*3]) == NULL) {
    printf("Error allocating memory !\n");
    return false;
  }
  

  // smooth the vectors on the surface
  //
  for(s=0; s < smSteps; s++) {  
    changed = false;
    if(smSteps > 1) {
      printf("Smoothing step %d.\r", s);
      fflush(stdout);
    }
    
    if(N >= 3) {
      // initialize buffer with first 3 slices of data (forces)
      for(k=0; k < 3; k++) {
	for (j = 0; j < M; j++) {
	  for (i = 0; i < L; i++) {
	    idx = k*slsz + j*L + i;
	    bidx = k*slsz + j*L + i;
	    
	    buffer[bidx].xd = gf_ws * force[idx].xd;
	    buffer[bidx].yd = gf_ws * force[idx].yd;
	    buffer[bidx].zd = gf_ws * force[idx].zd;
	  }
	}
      }
    }

    //
    // start smoothing
    for (k = 1; k < Nm1; k++) {

      // 
      // before processing each new slice, scroll the buffer
      if(k > 1) {
	
	// the first slice of the buffer is written back to slice k-2
	memcpy(force + ((k-2)*slsz), buffer, sizeof(ForceVector)*slsz);
	    
	// the last two slices are moved to the first two
	memcpy(buffer, buffer + slsz, sizeof(ForceVector)*slsz);
	memcpy(buffer + slsz, buffer + (2 * slsz), sizeof(ForceVector)*slsz);

	// the last slice is loaded from force from slice k+1
	memcpy(buffer + (2 * slsz), force + ((k + 1)*slsz), 
	       sizeof(ForceVector)*slsz);
	// multiply forces in this slice with self-weight
	bk = 2;
	for(idx=0; idx < L*M; idx++) {
	  bidx = bk*slsz + idx;
	  
	  buffer[bidx].xd = gf_ws * buffer[bidx].xd;
	  buffer[bidx].yd = gf_ws * buffer[bidx].yd;
	  buffer[bidx].zd = gf_ws * buffer[bidx].zd;
	}
      }

      //
      // processing the slice
      // values are read from the force array which has old values
      // and are written into the buffer which is then flushed when moving to 
      // the next slice
      for (j = 1; j < Mm1; j++) {
	for (i = 1; i < Lm1; i++) {  
	  idx = k*slsz + j*L + i;
	  bidx = slsz + j*L + i;
	  
	  if(f[idx] == SURF) {
	    // look at the neighbors and average the forces if not 0
	    //
	    
	    // save previous value
	    oxd = force[idx].xd;
	    oyd = force[idx].yd;
	    ozd = force[idx].zd;

	    // face neighbors
	    for(p=0; p < 6; p++) {
	      iidx = idx + ng[p];		// index of neighbor
	      // take only neighbors that are SURF
	      if(f[iidx] == SURF) {
		buffer[bidx].xd = buffer[bidx].xd + (gf_wf * force[iidx].xd);
		buffer[bidx].yd = buffer[bidx].yd + (gf_wf * force[iidx].yd);
		buffer[bidx].zd = buffer[bidx].zd + (gf_wf * force[iidx].zd);
	      }
	    }

	    // edge
	    for(p=6; p < 18; p++) {
	      iidx = idx + ng[p];		// index of neighbor
	      // take only neighbors that are SURF
	      if(f[iidx] == SURF) {
		buffer[bidx].xd = buffer[bidx].xd + (gf_we * force[iidx].xd);
		buffer[bidx].yd = buffer[bidx].yd + (gf_we * force[iidx].yd);
		buffer[bidx].zd = buffer[bidx].zd + (gf_we * force[iidx].zd);
	      }
	    }
	    
	    // vertex
	    for(p=18; p < 26; p++) {
	      iidx = idx + ng[p];		// index of neighbor
	      // take only neighbors that are SURF
	      if(f[iidx] == SURF) {
		buffer[bidx].xd = buffer[bidx].xd + (gf_wv * force[iidx].xd);
		buffer[bidx].yd = buffer[bidx].yd + (gf_wv * force[iidx].yd);
		buffer[bidx].zd = buffer[bidx].zd + (gf_wv * force[iidx].zd);
	      }
	    }

	    //
	    // normalize forces in buffer before writing them out
	    //
	    Normalize(buffer + bidx);

	    //	    
	    // record change
	    //
	    if((buffer[bidx].xd != oxd) || (buffer[bidx].yd != oyd) ||
	       (buffer[bidx].zd != ozd))
	    {
	      changed = true;
	      // add amount of change, computed as dot product of the original 
	      // and the new vector
	      change = 1 - 
		fabs((buffer[bidx].xd * oxd) + (buffer[bidx].yd * oyd) +
		     (buffer[bidx].zd * ozd));
	      
	      (*totalChange) = (*totalChange) + change;
	      if(change > (*maxChange)) {
		(*maxChange) = change;
	      }
	      else {
		if(change < (*minChange)) {
		  (*minChange) = change;
		}
	      }
	      (*nrChanged) = (*nrChanged) + 1;
	    }
	    else {
	      // printf("force unchanged !\n");
	    }
	  }
	}
      }
    }
    
    // the first and second slices of the buffer are written back to slices
    // N-3 and N-2
    if((N-3) >= 0) {
      memcpy(force + ((N-3) * slsz), buffer, sizeof(ForceVector)*slsz);
      memcpy(force + ((N-2) * slsz), buffer + slsz, sizeof(ForceVector)*slsz);
    }

    if(!changed) {
      // nothing changed in this step - stop
      printf("Smoothing step %d. - nothing changed.\n", s);
      fflush(stdout);
      break;
    }
  }
  
  delete [] buffer;
  return true;
}
			  



inline bool CheckNormals(unsigned char *fc, int L, int M, int N, 
		  ForceVector *force, int ng[26]) {
  int i, j, k, ni, nidx, idx;
  int Lm1, Mm1, Nm1, slsz;
  double r;
  
  Lm1 = L - 1;
  Mm1 = M - 1;
  Nm1 = N - 1;

  slsz = L*M;

  // assumes surface vectors are normalized 
  //printf("Checking surface vectors...");
  //fflush(stdout);
  for (k = 1; k < Nm1; k++) {
    for (j = 1; j < Mm1; j++) {
      for (i = 1; i < Lm1; i++) {
	
	idx = k*slsz + j*L + i;

	if(fc[idx] == SURF) {
	  // ignore forces that are 0.00
	  if((!IS_ZERO(force[idx].xd)) || (!IS_ZERO(force[idx].yd)) ||
	     (!IS_ZERO(force[idx].zd)))
	  {
	    for(ni = 0; ni < 26; ni++) {
	      nidx = idx + ng[ni];
	      if(fc[nidx] == SURF) {
		// ignore neighbors that are 0.00
		if((!IS_ZERO(force[nidx].xd)) || (!IS_ZERO(force[nidx].yd)) ||
		   (!IS_ZERO(force[nidx].zd)))
		{
		  // compute dot product
		  r = force[idx].xd * force[nidx].xd + 
		    force[idx].yd * force[nidx].yd + 
		    force[idx].zd * force[nidx].zd;
		  if(EQUAL(r, 1.00)) {
		    // these 2 vectors are the same.
		    return false;
		  }
		}
	      }
	    }
	  }
	}
      }
    }
  }
  //printf("done.\n");
  return true;

}





bool CriticalPointsAreClose(
  unsigned char *flags, 
  ForceVector *ForceField,
  unsigned char *cpMap,
  int *ng,
  int L, int M, int N) 
{
  int i, j, k, ii, slsz;
  int idx, cpIdx;
  int inds[8];
  CriticalPoint critPt;
  bool skipThisPoint;

  slsz = L*M;
  
  // clear critical point map
  memset(cpMap, 0, sizeof(unsigned char)*L*M*N);
     
  //
  for (k = 1; k < N-1; k++) {
    printf("\tProcessing plane %d out of %d\r", k, N-1);
    fflush(stdout);
    
    for (j = 1; j < M-1; j++) {
      for (i = 1; i < L-1; i++) {
	
	idx = k*slsz + j*L +i;
	
	// these are the neighbors that will be touched by the interpolation
	inds[0] = idx;
	inds[1] = idx + 1;
	inds[2] = idx + L;
	inds[3] = idx + L + 1;
	inds[4] = idx + slsz;
	inds[5] = idx + slsz + 1;
	inds[6] = idx + slsz + L;
	inds[7] = idx + slsz + L + 1;

	skipThisPoint = false;

	// we should skip the cells that have a vertex outside the object
	// if we don't, we will end up finding critical points outside the 
	// object !!
	for(ii=0; ii < 8; ii++) {
	  if(flags[inds[ii]] == EXTERIOR) {
	    skipThisPoint = true;
	    break;
	  }
	}
	
	if(skipThisPoint)  {
	  continue;
	}

	//
	// not skipped
	//
	
	if(FindCriticalPointInIntCell(i, j, k, inds, L, M, N, 
				      ForceField, flags, false, 
				      10, true, 
				      &(critPt.position))) 
	{
	  // classify critical point
	  if(ClassifyCriticalPoint(critPt, ForceField, flags, L, M, N)) {
	    if(critPt.type == CPT_ATTRACTING_NODE) {
	      cpIdx = (int(critPt.position.z)*L*M) + 
		(int(critPt.position.y)*L) + int(critPt.position.x);

	      // mark the critical point
	      cpMap[cpIdx] = 1;
	      
	      // check the neighbors 
	      for(ii=0; ii < 26; ii++) {
		if(cpMap[cpIdx + ng[ii]] != 0) {
		  // it has a critical point as a neighbor
		  return true;
		}
	      }
	      // not close
	    }
	  }
	}
      } // for i
    } // for j
  } // for k

  return false;
}




bool SetNeighborsArray(int ng[26], int L, int M, int N) {
  int slsz = L*M;		// slice size
  
  // face neighbors
  ng[0]	= + slsz + 0 + 0;
  ng[1]	= - slsz + 0 + 0;
  ng[2]	= +    0 + L + 0;
  ng[3]	= +    0 - L + 0;
  ng[4]	= +    0 + 0 + 1;
  ng[5]	= +    0 + 0 - 1;
  // e-neighbors
  ng[6]	= + slsz + L + 0;
  ng[7]	= + slsz - L + 0;
  ng[8]	= - slsz + L + 0;
  ng[9]	= - slsz - L + 0;
  ng[10]	= + slsz + 0 + 1;
  ng[11]	= + slsz + 0 - 1;
  ng[12]	= - slsz + 0 + 1;
  ng[13]	= - slsz + 0 - 1;
  ng[14]	= +    0 + L + 1;
  ng[15]	= +    0 + L - 1;
  ng[16]	= +    0 - L + 1;
  ng[17]	= +    0 - L - 1;
  // v-neighbors
  ng[18]	= - slsz - L - 1;
  ng[19]	= - slsz - L + 1;
  ng[20]	= - slsz + L - 1;
  ng[21]	= - slsz + L + 1;
  ng[22]	= + slsz - L - 1;
  ng[23]	= + slsz - L + 1;
  ng[24]	= + slsz + L - 1;
  ng[25]	= + slsz + L + 1;

  return true;
}




bool ComputeSurfaceNormals(unsigned char *fc, int L, int M, int N, 
			   ForceVector *force)
{
  int i, j, k, idx, slsz;
  int Lm1, Mm1, Nm1;
  bool hasNeighbors;

  Lm1 = L - 1;
  Mm1 = M - 1;
  Nm1 = N - 1;
  slsz = L*M;		// slice size

  // compute gradient of surface voxels
  for (k = 1; k < Nm1; k++) {
    for (j = 1; j < Mm1; j++) {
      for (i = 1; i < Lm1; i++) {
	
	idx = k*slsz + j*L + i;
	
	if(fc[idx] == SURF) {
	  // surface voxel
	  hasNeighbors = false;
	  
	  force[idx].xd = 0;
	  if(fc[idx + 1] != EXTERIOR)  {
	    hasNeighbors = true;
	    force[idx].xd = force[idx].xd + 1;
	  }
	  if(fc[idx - 1] != EXTERIOR)  {
	    hasNeighbors = true;
	    force[idx].xd = force[idx].xd - 1;
	  }
	  
	  force[idx].yd = 0;
	  if(fc[idx + L] != EXTERIOR)  {
	    hasNeighbors = true;
	    force[idx].yd = force[idx].yd + 1;
	  }
	  if(fc[idx - L] != EXTERIOR)  {
	    hasNeighbors = true;
	    force[idx].yd = force[idx].yd - 1;
	  }

	  force[idx].zd = 0;
	  if(fc[idx + slsz] != EXTERIOR)  {
	    hasNeighbors = true;
	    force[idx].zd = force[idx].zd + 1;
	  }
	  if(fc[idx - slsz] != EXTERIOR)  {
	    hasNeighbors = true;
	    force[idx].zd = force[idx].zd - 1;
	  }
	  
	  if(!hasNeighbors) {
	    printf("WARNING: Surface voxel has no neighbors. Force unchanged.\n");	    
	  }
	}
      }
    }
  }
  
  return true;
}




bool ComputeSubSurfaceVectors(ForceVector *force, unsigned char *fc, 
			      int L, int M, int N, int ng[26])
{
  bool done;
  int cnt;
  int i, j, k, s, idx, iidx, slsz;
  
  slsz = L*M;
  
  done = false;
  cnt = 1;
  while(!done) {    
    cnt++;
    printf("layer %d ...", cnt);
    fflush(stdout);
    
    done = true;
    // flag sub-surface voxels
    // look at the INTERIOR voxels. 
    // If an INTERIOR voxel has an SURF neighbor,
    //    then it is a SUB-SURFace voxel
    // Look only at the face neighbors.
    for (k = 1; k < (N-1); k++) {
      for (j = 1; j < (M-1); j++) {
	for (i = 1; i < (L-1); i++) {
	  idx = k*slsz + j*L + i;
	  
	  if(fc[idx] == INTERIOR) {
	    for(s=0; s < 6; s++) {
	      iidx = idx + ng[s];
	      
	      if(fc[iidx] == SURF) {
		fc[idx] = SUB_SURF;
		done = false;
		break;
	      }
	    }
	  }
	}
      }
    }
    
    if(done) {
      // there are no sub-surface voxels - nothing more to do
      printf("no more sub-surface voxels. - done\n");
      break;
    }

    // compute the force at sub-surface voxels by averaging the forces
    // at the surface neighbors
    for (k = 1; k < (N-1); k++) {
      for (j = 1; j < (M-1); j++) {
	for (i = 1; i < (L-1); i++) {
	  
	  idx = k*slsz + j*L + i;
	  
	  if(fc[idx] == SUB_SURF) {
	    // look at the neighbors and average the forces if not 0
	    //
	    force[idx].xd = 0;
	    force[idx].yd = 0;
	    force[idx].zd = 0;

	    // face
	    for(s=0; s < 6; s++) {
	      iidx = idx + ng[s];		// index of neighbor
	      // take only neighbors that are SURF
	      if(fc[iidx] == SURF) {
		force[idx].xd = force[idx].xd + (gf_wf * force[iidx].xd);
		force[idx].yd = force[idx].yd + (gf_wf * force[iidx].yd);
		force[idx].zd = force[idx].zd + (gf_wf * force[iidx].zd);
	      }
	    }
	    
	    // edge
	    for(s=6; s < 18; s++) {
	      iidx = idx + ng[s];		// index of neighbor
	      // take only neighbors that are SURF
	      if(fc[iidx] == SURF) {
		force[idx].xd = force[idx].xd + (gf_we * force[iidx].xd);
		force[idx].yd = force[idx].yd + (gf_we * force[iidx].yd);
		force[idx].zd = force[idx].zd + (gf_we * force[iidx].zd);
	      }
	    }
	    
	    for(s=18; s < 26; s++) {
	      iidx = idx + ng[s];		// index of neighbor
	      // take only neighbors that are SURF
	      if(fc[iidx] == SURF) {
		force[idx].xd = force[idx].xd + (gf_wv * force[iidx].xd);
		force[idx].yd = force[idx].yd + (gf_wv * force[iidx].yd);
		force[idx].zd = force[idx].zd + (gf_wv * force[iidx].zd);
	      }
	    }
	    
	    // normalize vector
	    Normalize(force + idx);
	    
	  }
	}
      }
    }

    
    // remove surface voxels
    for(idx = 0; idx < (L*M*N); idx++) {
      if(fc[idx] == SURF) {
	fc[idx] = EXTERIOR;
      }
    }
    
    // re-flag interior, boundary and exterior voxels
    FlagVolume(fc, L, M, N);
    
    printf("done.\r");
    fflush(stdout);
  }
  
  return true;
}
