#include "../src/grm.h"
#include "mex.h"
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  auto geno = reinterpret_cast<uint8_t *>(mxGetData(prhs[0]));
  uint32_t V = static_cast<uint32_t>(mxGetN(prhs[0]));
  auto af = mxGetPr(prhs[1]);
  uint32_t N = static_cast<uint32_t>(mxGetScalar(prhs[2]));
  uint32_t num_threads = static_cast<uint32_t>(mxGetScalar(prhs[3]));
  plhs[0] = mxCreateDoubleMatrix(N, N, mxREAL);
  auto grm = mxGetPr(plhs[0]);
  CalcGrmMatrix(geno, af, N, V, grm, num_threads);
}
