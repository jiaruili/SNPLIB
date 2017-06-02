#ifndef SNPLIB_BASE_H
#define SNPLIB_BASE_H

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

void CalcAlleleFrequency(const uint8_t *geno, uint32_t num_indivs,
                         uint32_t num_snps, double *af, uint32_t num_threads);
void CalcCallrate(const uint8_t *geno, uint32_t num_indivs, uint32_t num_snps,
                  double *cr, uint32_t num_threads);
void SubsetGeno(const uint8_t *geno_src, uint32_t num_indivs,
                uint32_t *snp_index, uint32_t num_snps, uint8_t *geno_dest,
                uint32_t num_threads);

#ifdef __cplusplus
}
#endif
#endif /* end of include guard: SNPLIB_BASE_H */
