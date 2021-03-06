#ifndef SNPLIB_GRM_H
#define SNPLIB_GRM_H

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

void CalcGrmMatrix(const uint8_t *geno, const double *af, uint32_t num_indivs,
                   uint32_t num_snps, double *matrix, uint32_t num_threads);
#ifdef __cplusplus
}
#endif
#endif /* end of include guard: SNPLIB_GRM_H */
