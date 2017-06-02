#include "../src/spectral_ibs.h"
#include "mex.h"
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  auto geno_dest = reinterpret_cast<uint8_t *>(mxGetData(prhs[0]));
  auto geno_src = reinterpret_cast<uint8_t *>(mxGetData(prhs[1]));
  uint32_t V = static_cast<uint32_t>(mxGetN(prhs[0]));
  double *loadings = mxGetPr(prhs[2]);
  uint32_t K = static_cast<uint32_t>(mxGetN(prhs[2]));
  uint32_t N_src = static_cast<uint32_t>(mxGetScalar(prhs[3]));
  uint32_t N_dest = static_cast<uint32_t>(mxGetScalar(prhs[4]));
  uint32_t num_threads = static_cast<uint32_t>(mxGetScalar(prhs[5]));
  plhs[0] = mxCreateDoubleMatrix(N_dest, K, mxREAL);
  auto scores = mxGetPr(plhs[0]);
  ProjectSpectralIBSSpace(geno_src, geno_dest, loadings, N_src, N_dest, V, K,
                          scores, num_threads);
}
