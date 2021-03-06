//
// compute distance transform
//
bool GetDT(
  unsigned char *cf,    // [in]  input volume
  int L, int M, int N,  // [in]  volume size
  float *f              // [out] distance field
);


///////////////////////////////////////////////////////////////////////////////
// reads distance field data from a given file
///////////////////////////////////////////////////////////////////////////////
bool ReadDistanceField(char *filename, int L, int M, int N, float **df);

///////////////////////////////////////////////////////////////////////////////
// write distance field data from a given file
///////////////////////////////////////////////////////////////////////////////
bool SaveDistanceField(char *filename, int L, int M, int N, float *df);
