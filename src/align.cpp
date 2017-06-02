#include "align.h"

namespace {
const uint32_t snps_byte = 4;
const uint32_t snps_word = 32;
void FlipThread(uint8_t *geno, uint32_t bytes_snp, uint32_t words_snp,
                uint32_t num_snps, uint32_t *snps_list) {
  uint64_t *geno64 = new uint64_t[words_snp];
  memset(geno64, 0x55u, sizeof(uint64_t) * words_snp);
  for (uint32_t i = 0; i < num_snps; ++i) {
    auto g8 = geno + snps_list[i] * bytes_snp;
    memcpy(geno64, g8, sizeof(uint8_t) * bytes_snp);
    for (uint32_t j = 0; j < words_snp; ++j) {
      uint64_t tmp_geno = geno64[j] ^ 0x5555555555555555ull;
      uint64_t tmp_mask_1 =
          (tmp_geno | (tmp_geno >> 1)) & 0x5555555555555555ull;
      tmp_mask_1 *= 3;
      tmp_geno = (~geno64[j]) ^ 0x5555555555555555ull;
      uint64_t tmp_mask_2 =
          (tmp_geno | (tmp_geno >> 1)) & 0x5555555555555555ull;
      tmp_mask_2 *= 3;
      uint64_t mask = tmp_mask_1 | tmp_mask_2;
      geno64[j] = (geno64[j] & mask) + ~(geno64[j] | mask);
    }
    memcpy(g8, geno64, sizeof(uint8_t) * bytes_snp);
  }
}
} // namespace

void FlipGeno(uint8_t *geno, uint32_t num_indivs, uint32_t num_snps,
              uint32_t *snps_list, uint32_t num_threads) {
  uint32_t bytes_snp =
      num_indivs / snps_byte + (num_indivs % snps_byte ? 1 : 0);
  uint32_t words_snp =
      num_indivs / snps_word + (num_indivs % snps_word ? 1 : 0);
  std::vector<std::thread> workers(num_threads);
  uint32_t snps_job = num_snps / num_threads + 1;
  uint32_t snps_left = num_snps % num_threads;
  for (uint32_t i = 0; i < snps_left; ++i) {
    workers[i] = std::thread(FlipThread, geno, bytes_snp, words_snp, snps_job,
                             snps_list);
    snps_list += snps_job;
  }
  --snps_job;
  for (uint32_t i = snps_left; i < num_threads; ++i) {
    workers[i] = std::thread(FlipThread, geno, bytes_snp, words_snp, snps_job,
                             snps_list);
    snps_list += snps_job;
  }
  for (auto &&iter : workers) {
    iter.join();
  }
}
