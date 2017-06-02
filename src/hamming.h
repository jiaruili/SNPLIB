#ifndef SNPLIB_HAMMING_H
#define SNPLIB_HAMMING_H

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

void CalcHammingMatrix(const uint8_t *geno, uint32_t num_indivs,
                       uint32_t num_snps, double *matrix, uint32_t num_threads);
#ifdef __cplusplus
}
#endif
#endif /* end of include guard: SNPLIB_HAMMING_H */
