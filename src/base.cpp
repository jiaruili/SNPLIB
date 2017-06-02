#include "base.h"

namespace {
const uint32_t snps_byte = 4;
const uint32_t snps_word = 32;
void CalcAfThread(const uint8_t *geno, uint32_t bytes_snp, uint32_t words_snp,
                  uint32_t correct, uint32_t num_snps, double *af) {
  uint64_t *geno64 = new uint64_t[words_snp];
  memset(geno64, 0x55u, sizeof(uint64_t) * words_snp);
  for (uint32_t i = 0; i < num_snps; ++i) {
    auto g8 = geno + i * bytes_snp;
    memcpy(geno64, g8, sizeof(uint8_t) * bytes_snp);
    uint64_t alleles = 0;
    uint64_t nonmissings = 0;
    for (uint32_t j = 0; j < words_snp; ++j) {
      uint64_t tmp_geno = geno64[j] ^ 0x5555555555555555ull;
      uint64_t tmp_mask = (tmp_geno | (tmp_geno >> 1)) & 0x5555555555555555ull;
      tmp_mask *= 3;
      nonmissings += _mm_popcnt_u64(tmp_mask);
      tmp_geno = tmp_mask & geno64[j];
      alleles += _mm_popcnt_u64(tmp_geno);
    }
    nonmissings -= correct;
    af[i] = static_cast<double>(alleles) / nonmissings;
  }
}

void CalcCrThread(const uint8_t *geno, uint32_t bytes_snp, uint32_t words_snp,
                  uint32_t correct, uint32_t num_snps, double *cr) {
  uint64_t *geno64 = new uint64_t[words_snp];
  memset(geno64, 0x55u, sizeof(uint64_t) * words_snp);
  for (uint32_t i = 0; i < num_snps; ++i) {
    auto g8 = geno + i * bytes_snp;
    memcpy(geno64, g8, sizeof(uint8_t) * bytes_snp);
    uint64_t nonmissings = 0;
    for (uint32_t j = 0; j < words_snp; ++j) {
      uint64_t tmp_geno = geno64[j] ^ 0x5555555555555555ull;
      uint64_t tmp_mask = (tmp_geno | (tmp_geno >> 1)) & 0x5555555555555555ull;
      tmp_mask *= 3;
      nonmissings += _mm_popcnt_u64(tmp_mask);
    }
    nonmissings -= correct;
    cr[i] = static_cast<double>(nonmissings);
  }
}
} // namespace

void CalcAlleleFrequency(const uint8_t *geno, uint32_t num_indivs,
                         uint32_t num_snps, double *af, uint32_t num_threads) {
  uint32_t bytes_snp =
      num_indivs / snps_byte + (num_indivs % snps_byte ? 1 : 0);
  uint32_t words_snp =
      num_indivs / snps_word + (num_indivs % snps_word ? 1 : 0);
  uint32_t indivs_left = num_indivs % snps_byte;
  uint32_t correct = indivs_left ? (snps_byte - indivs_left) : 0;
  correct *= 2;
  std::vector<std::thread> workers(num_threads);
  uint32_t snps_job = num_snps / num_threads + 1;
  uint32_t snps_left = num_snps % num_threads;
  for (uint32_t i = 0; i < snps_left; ++i) {
    workers[i] = std::thread(CalcAfThread, geno, bytes_snp, words_snp, correct,
                             snps_job, af);
    geno += snps_job * bytes_snp;
    af += snps_job;
  }
  --snps_job;
  for (uint32_t i = snps_left; i < num_threads; ++i) {
    workers[i] = std::thread(CalcAfThread, geno, bytes_snp, words_snp, correct,
                             snps_job, af);
    geno += snps_job * bytes_snp;
    af += snps_job;
  }
  for (auto &&iter : workers) {
    iter.join();
  }
}

void CalcCallrate(const uint8_t *geno, uint32_t num_indivs, uint32_t num_snps,
                  double *cr, uint32_t num_threads) {
  uint32_t bytes_snp =
      num_indivs / snps_byte + (num_indivs % snps_byte ? 1 : 0);
  uint32_t words_snp =
      num_indivs / snps_word + (num_indivs % snps_word ? 1 : 0);
  uint32_t indivs_left = num_indivs % snps_byte;
  uint32_t correct = indivs_left ? (snps_byte - indivs_left) : 0;
  correct *= 2;
  std::vector<std::thread> workers(num_threads);
  uint32_t snps_job = num_snps / num_threads + 1;
  uint32_t snps_left = num_snps % num_threads;
  for (uint32_t i = 0; i < snps_left; ++i) {
    workers[i] = std::thread(CalcCrThread, geno, bytes_snp, words_snp, correct,
                             snps_job, cr);
    geno += snps_job * bytes_snp;
    cr += snps_job;
  }
  --snps_job;
  for (uint32_t i = snps_left; i < num_threads; ++i) {
    workers[i] = std::thread(CalcCrThread, geno, bytes_snp, words_snp, correct,
                             snps_job, cr);
    geno += snps_job * bytes_snp;
    cr += snps_job;
  }
  for (auto &&iter : workers) {
    iter.join();
  }
  for (uint32_t i = 0; i < num_snps; ++i) {
    cr[i] /= num_indivs;
  }
}

void SubsetGeno(const uint8_t *geno_src, uint32_t num_indivs,
                uint32_t *snp_index, uint32_t num_snps, uint8_t *geno_dest,
                uint32_t num_threads) {}