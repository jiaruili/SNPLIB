#ifndef SNPLIB_IBS_H
#define SNPLIB_IBS_H
#include "barrier.h"
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

void CalcIbsMatrix(const uint8_t *geno, uint32_t num_indivs, uint32_t num_snps,
                   double *matrix, uint32_t num_threads);
void CalcIbsConnection(const uint8_t *geno_src, uint32_t num_indivs_src,
                       const uint8_t *geno_dest, uint32_t num_indivs_dest,
                       uint32_t num_snps, double *connection,
                       uint32_t num_threads);
#ifdef __cplusplus
}
#endif
#endif
