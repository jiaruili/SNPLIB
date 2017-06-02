#include "../src/base.h"
#include "mex.h"
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  auto geno = reinterpret_cast<uint8_t *>(mxGetData(prhs[0]));
  uint32_t V = static_cast<uint32_t>(mxGetN(prhs[0]));
  uint32_t N = static_cast<uint32_t>(mxGetScalar(prhs[1]));
  uint32_t num_threads = static_cast<uint32_t>(mxGetScalar(prhs[2]));
  plhs[0] = mxCreateDoubleMatrix(V, 1, mxREAL);
  auto AF = mxGetPr(plhs[0]);
  CalcAlleleFrequency(geno, N, V, AF, num_threads);
}
