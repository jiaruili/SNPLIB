#ifndef SNPLIB_PCA_H
#define SNPLIB_PCA_H

#include <cmath>
#include <cstring>
#include <thread>
#include <vector>
#ifdef _MSC_VER
#include <intrin.h>
#else
#include <x86intrin.h>
#endif

#ifdef __cplusplus
extern "C" {
#endif

void UProductPca(const uint8_t *geno, const double *af, uint32_t num_indivs,
                 uint32_t num_snps, uint32_t num_components, const double *B,
                 double *C, uint32_t num_threads);

void TProductPca(const uint8_t *geno, const double *af, uint32_t num_indivs,
                 uint32_t num_snps, uint32_t num_components, const double *B,
                 double *C, uint32_t num_threads);

#ifdef __cplusplus
}
#endif

#endif /* end of include guard: SNPLIB_PCA_H */
