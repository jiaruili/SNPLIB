#ifndef SNPLIB_ALIGN_H
#define SNPLIB_ALIGN_H

#include <cstring>
#include <thread>
#include <vector>
#ifdef __cplusplus
extern "C" {
#endif

void FlipGeno(uint8_t *geno, uint32_t num_indivs, uint32_t num_snps,
              uint32_t *snps_list, uint32_t num_threads);

#ifdef __cplusplus
}
#endif
#endif /* end of include guard: SNPLIB_ALIGN_H */
