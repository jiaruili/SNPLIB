#include "pca.h"

namespace {
const uint32_t snps_byte = 4;
const uint32_t snps_row = 384;
const uint32_t indivs_row = 384;
const uint32_t MC = snps_row;
const uint32_t KC = indivs_row;
const uint32_t NC = 256;
const uint32_t MR = 4;
const uint32_t NR = 4;
//
//  Packing complete panels from B (i.e. without padding)
//
static void pack_kxNR(uint32_t num_counts, const double *B, int ldb,
                      double *buffer) {
  for (uint32_t i = 0; i < num_counts; ++i) {
    for (uint32_t j = 0; j < NR; ++j) {
      buffer[j] = B[j * ldb];
    }
    buffer += NR;
    B++;
  }
}
//
//  Packing panels from B with padding if required
//
static void pack_B(int num_counts, int num_cols, const double *B, int ldb,
                   double *buffer) {
  uint32_t np = num_cols / NR;
  uint32_t _nr = num_cols % NR;
  for (uint32_t j = 0; j < np; ++j) {
    pack_kxNR(num_counts, B, ldb, buffer);
    buffer += num_counts * NR;
    B += NR * ldb;
  }
  if (_nr > 0) {
    for (uint32_t i = 0; i < num_counts; ++i) {
      for (uint32_t j = 0; j < _nr; ++j) {
        buffer[j] = B[j * ldb];
      }
      for (uint32_t j = _nr; j < NR; ++j) {
        buffer[j] = 0.0;
      }
      buffer += NR;
      B++;
    }
  }
}
//
//  Micro kernel for multiplying panels from A and B.
//
static void dgemm_micro_kernel(uint32_t num_counts, const double *A,
                               const double *B, double beta, double *C,
                               uint32_t ldc) {
  double AB[MR * NR] __attribute__((aligned(16)));
  //
  //  Compute AB = A*B
  //
  __m128d ab_00_11, ab_20_31;
  __m128d ab_01_10, ab_21_30;
  __m128d ab_02_13, ab_22_33;
  __m128d ab_03_12, ab_23_32;

  __m128d tmp0, tmp1, tmp2, tmp3;
  __m128d tmp4, tmp5, tmp6, tmp7;

  tmp0 = _mm_load_pd(A);     // (1)
  tmp1 = _mm_load_pd(A + 2); // (2)
  tmp2 = _mm_load_pd(B);     // (3)

  ab_00_11 = _mm_setzero_pd();
  ab_20_31 = _mm_setzero_pd();
  ab_01_10 = _mm_setzero_pd();
  ab_21_30 = _mm_setzero_pd();
  ab_02_13 = _mm_setzero_pd();
  ab_22_33 = _mm_setzero_pd();
  ab_03_12 = _mm_setzero_pd();
  ab_23_32 = _mm_setzero_pd();

  for (uint32_t l = 0; l < num_counts; ++l) {
    tmp3 = _mm_load_pd(B + 2);
    tmp4 = _mm_shuffle_pd(tmp2, tmp2, _MM_SHUFFLE2(0, 1));
    tmp5 = _mm_shuffle_pd(tmp3, tmp3, _MM_SHUFFLE2(0, 1));
    tmp6 = tmp2;
    tmp2 = _mm_mul_pd(tmp2, tmp0);
    tmp6 = _mm_mul_pd(tmp6, tmp1);
    ab_00_11 = _mm_add_pd(ab_00_11, tmp2);
    ab_20_31 = _mm_add_pd(ab_20_31, tmp6);
    tmp7 = tmp4;
    tmp4 = _mm_mul_pd(tmp4, tmp0);
    tmp7 = _mm_mul_pd(tmp7, tmp1);
    ab_01_10 = _mm_add_pd(ab_01_10, tmp4);
    ab_21_30 = _mm_add_pd(ab_21_30, tmp7);
    tmp2 = _mm_load_pd(B + 4); // (6)
    tmp6 = tmp3;
    tmp3 = _mm_mul_pd(tmp3, tmp0);
    tmp6 = _mm_mul_pd(tmp6, tmp1);
    ab_02_13 = _mm_add_pd(ab_02_13, tmp3);
    ab_22_33 = _mm_add_pd(ab_22_33, tmp6);
    tmp7 = tmp5;
    tmp5 = _mm_mul_pd(tmp5, tmp0);
    tmp0 = _mm_load_pd(A + 4); // (4)
    tmp7 = _mm_mul_pd(tmp7, tmp1);
    tmp1 = _mm_load_pd(A + 6); // (5)
    ab_03_12 = _mm_add_pd(ab_03_12, tmp5);
    ab_23_32 = _mm_add_pd(ab_23_32, tmp7);
    A += 4;
    B += 4;
  }
  _mm_storel_pd(&AB[0 + 0 * 4], ab_00_11);
  _mm_storeh_pd(&AB[1 + 0 * 4], ab_01_10);
  _mm_storel_pd(&AB[2 + 0 * 4], ab_20_31);
  _mm_storeh_pd(&AB[3 + 0 * 4], ab_21_30);
  _mm_storel_pd(&AB[0 + 1 * 4], ab_01_10);
  _mm_storeh_pd(&AB[1 + 1 * 4], ab_00_11);
  _mm_storel_pd(&AB[2 + 1 * 4], ab_21_30);
  _mm_storeh_pd(&AB[3 + 1 * 4], ab_20_31);

  _mm_storel_pd(&AB[0 + 2 * 4], ab_02_13);
  _mm_storeh_pd(&AB[1 + 2 * 4], ab_03_12);
  _mm_storel_pd(&AB[2 + 2 * 4], ab_22_33);
  _mm_storeh_pd(&AB[3 + 2 * 4], ab_23_32);
  _mm_storel_pd(&AB[0 + 3 * 4], ab_03_12);
  _mm_storeh_pd(&AB[1 + 3 * 4], ab_02_13);
  _mm_storel_pd(&AB[2 + 3 * 4], ab_23_32);
  _mm_storeh_pd(&AB[3 + 3 * 4], ab_22_33);
  //
  //  Update C <- beta*C
  //
  if (beta == 0.0) {
    for (uint32_t j = 0; j < NR; ++j) {
      for (uint32_t i = 0; i < MR; ++i) {
        C[i + j * ldc] = 0.0;
      }
    }
  } else if (beta != 1.0) {
    for (uint32_t j = 0; j < NR; ++j) {
      for (uint32_t i = 0; i < MR; ++i) {
        C[i + j * ldc] *= beta;
      }
    }
  }
  //
  //  Update C <- C + alpha*AB (note: the case alpha==0.0 was already treated in
  //                                  the above layer dgemm_nn)
  //
  for (uint32_t j = 0; j < NR; ++j) {
    for (uint32_t i = 0; i < MR; ++i) {
      C[i + j * ldc] += AB[i + j * MR];
    }
  }
}
//
//  Compute Y += alpha*X
//
static void dgeaxpy(uint32_t m, uint32_t n, double alpha, const double *X,
                    uint32_t ldx, double *Y, uint32_t ldy) {
  int i, j;

  if (alpha != 1.0) {
    for (j = 0; j < n; ++j) {
      for (i = 0; i < m; ++i) {
        Y[i + j * ldy] += alpha * X[i + j * ldx];
      }
    }
  } else {
    for (j = 0; j < n; ++j) {
      for (i = 0; i < m; ++i) {
        Y[i + j * ldy] += X[i + j * ldx];
      }
    }
  }
}
//
//  Compute X *= alpha
//
static void dgescal(uint32_t m, uint32_t n, double alpha, double *X, int ldx) {
  int i, j;
  if (alpha != 0.0) {
    for (j = 0; j < n; ++j) {
      for (i = 0; i < m; ++i) {
        X[i + j * ldx] *= alpha;
      }
    }
  } else {
    for (j = 0; j < n; ++j) {
      for (i = 0; i < m; ++i) {
        X[i + j * ldx] = 0.0;
      }
    }
  }
}
//
//  Macro Kernel for the multiplication of blocks of A and B.  We assume that
//  these blocks were previously packed to buffers _A and _B.
//
static void dgemm_macro_kernel(uint32_t num_rows, uint32_t num_cols,
                               uint32_t num_counts, double beta,
                               double *buffer_a, double *buffer_b, double *C,
                               uint32_t ldc) {
  uint32_t mp = (num_rows + MR - 1) / MR;
  uint32_t np = (num_cols + NR - 1) / NR;
  uint32_t _mr = num_rows % MR;
  uint32_t _nr = num_cols % NR;
  int mr, nr;
  int i, j;
  double buffer_c[MR * NR];
  for (j = 0; j < np; ++j) {
    nr = (j != np - 1 || _nr == 0) ? NR : _nr;
    for (i = 0; i < mp; ++i) {
      mr = (i != mp - 1 || _mr == 0) ? MR : _mr;
      if (mr == MR && nr == NR) {
        dgemm_micro_kernel(num_counts, &buffer_a[i * num_counts * MR],
                           &buffer_b[j * num_counts * NR], beta,
                           &C[i * MR + j * NR * ldc], ldc);
      } else {
        dgemm_micro_kernel(num_counts, &buffer_a[i * num_counts * MR],
                           &buffer_b[j * num_counts * NR], 0.0, buffer_c, MR);
        dgescal(mr, nr, beta, &C[i * MR + j * NR * ldc], ldc);
        dgeaxpy(mr, nr, 1.0, buffer_c, MR, &C[i * MR + j * NR * ldc], ldc);
      }
    }
  }
}

void UpdateAf(const double af, double std_af[4]) {
  double std = sqrt(2.0 * af * (1 - af));
  std_af[0] = -2.0 * af / std;
  std_af[1] = 0.0;
  std_af[2] = (1.0 - 2.0 * af) / std;
  std_af[3] = (2.0 - 2.0 * af) / std;
}

void UnfoldGenoUFull(const uint8_t *geno, const double *std_af,
                     uint32_t num_indivs, uint32_t bytes_snp,
                     uint32_t full_bytes, double *geno_d) {
  uint32_t blocks_row = snps_row / MR;
  for (uint32_t i = 0; i < blocks_row; ++i) {
    auto tmp_g = geno + i * MR * bytes_snp;
    auto tmp_d = geno_d + i * MR * num_indivs;
    auto tmp_af = std_af + 4 * i * MR;
    for (uint32_t k = 0; k < MR; ++k) {
      auto g = tmp_g + k * bytes_snp;
      auto d = tmp_d + k;
      auto af = tmp_af + 4 * k;
      for (uint32_t l = 0; l < full_bytes; ++l) {
        uint8_t t = g[l];
        d[4 * l * MR] = af[t & 3u];
        t >>= 2;
        d[(4 * l + 1) * MR] = af[t & 3u];
        t >>= 2;
        d[(4 * l + 2) * MR] = af[t & 3u];
        t >>= 2;
        d[(4 * l + 3) * MR] = af[t & 3u];
      }
    }
  }
}

void UProdPcaThread(const uint8_t *geno, const double *af, uint32_t num_indivs,
                    uint32_t num_snps, uint32_t num_blocks,
                    uint32_t num_components, const double *B, double *C) {
  uint32_t full_bytes = num_indivs / snps_byte;
  uint32_t indivs_left = num_indivs % snps_byte;
  uint32_t bytes_snp = full_bytes + (indivs_left ? 1 : 0);
  uint32_t kb = num_indivs / KC;
  uint32_t _kc = num_indivs % KC;
  uint32_t nc = num_components;
  double beta;
  double *geno_d = new double[snps_row * KC];
  double *std_af = new double[snps_row * 4];
  double *buffer_b = new double[KC * NC];
  for (uint32_t l = 0; l < num_blocks; ++l) {
    auto tmp_g = geno + l * snps_row * bytes_snp;
    auto tmp_c = C + l * snps_row;
    auto tmp_af = af + snps_row;
    for (uint32_t i = 0; i < snps_row; ++i) {
      UpdateAf(tmp_af[i], std_af + 4 * i);
    }
    for (uint32_t i = 0; i < kb; ++i) {
      beta = i == 0 ? 0.0 : 1.0;
      UnfoldGenoUFull(tmp_g + i * KC / snps_byte, std_af, KC, bytes_snp,
                      full_bytes, geno_d);
      pack_B(KC, nc, B + i * KC, num_indivs, buffer_b);
      dgemm_macro_kernel(MC, nc, KC, beta, geno_d, buffer_b, tmp_c, num_snps);
    }
  }
  delete[] geno_d;
  delete[] buffer_b;
  delete[] std_af;
}
} // namespace

void UProductPca(const uint8_t *geno, const double *af, uint32_t num_indivs,
                 uint32_t num_snps, uint32_t num_components, const double *B,
                 double *C, uint32_t num_threads) {}

void TProductPca(const uint8_t *geno, const double *af, uint32_t num_indivs,
                 uint32_t num_snps, uint32_t num_components, const double *B,
                 double *C, uint32_t num_threads) {}
