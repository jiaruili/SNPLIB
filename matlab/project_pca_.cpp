#include "../src/pca.h"
#include "mex.h"
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  auto geno_dest = reinterpret_cast<uint8_t *>(mxGetData(prhs[0]));
  uint32_t V = static_cast<uint32_t>(mxGetN(prhs[0]));
  double *af = mxGetPr(prhs[1]);
  double *loadings = mxGetPr(prhs[2]);
  uint32_t K = static_cast<uint32_t>(mxGetN(prhs[2]));
  uint32_t N = static_cast<uint32_t>(mxGetScalar(prhs[3]));
  uint32_t num_threads = static_cast<uint32_t>(mxGetScalar(prhs[4]));
  plhs[0] = mxCreateDoubleMatrix(N, K, mxREAL);
  auto scores = mxGetPr(plhs[0]);
  ProjectPcaSpace(geno_dest, af, loadings, N, V, K, scores, num_threads);
}
