#include "hamming.h"

namespace {
const uint32_t snps_byte = 4;
const uint32_t snps_word = 20;
const uint32_t words_block = 32;
const uint64_t trans_mask = 0x0003000300030003ull; // 0003 = 0000000000000011
const uint64_t mask_mask = 0x1249124912491249ull;  // 1249 = 0001001001001001
const auto snps_block = snps_word * words_block;
const auto dwords_block = words_block / 2;

void ConvertGeno(const uint8_t *geno, uint32_t bytes_snp, uint64_t tmp_g64[5]) {
  uint64_t tmp_g8[snps_word];
  for (uint32_t i = 0; i < snps_word; ++i) {
    tmp_g8[i] = static_cast<uint64_t>(geno[i * bytes_snp]);
  }
  for (uint32_t i = 0; i < 5; ++i) {
    tmp_g64[i] = tmp_g8[i] + (tmp_g8[5 + i] << 16) + (tmp_g8[10 + i] << 32) +
                 (tmp_g8[15 + i] << 48);
  }
}

void ConvertLastGeno(const uint8_t *geno, uint32_t bytes_snp,
                     uint32_t snps_left, uint64_t tmp_g64[5]) {
  uint64_t tmp_g8[snps_word];
  for (uint32_t i = 0; i < snps_left; ++i) {
    tmp_g8[i] = static_cast<uint64_t>(geno[i * bytes_snp]);
  }
  for (uint32_t i = snps_left; i < snps_word; ++i) {
    tmp_g8[i] = 0x55ull;
  }
  for (uint32_t i = 0; i < 5; ++i) {
    tmp_g64[i] = tmp_g8[i] + (tmp_g8[5 + i] << 16) + (tmp_g8[10 + i] << 32) +
                 (tmp_g8[15 + i] << 48);
  }
}

uint64_t ConvertG64(const uint64_t tmp_g64[5]) {
  uint64_t geno64;
  geno64 = tmp_g64[4] & trans_mask;
  geno64 <<= 3;
  geno64 |= tmp_g64[3] & trans_mask;
  geno64 <<= 3;
  geno64 |= tmp_g64[2] & trans_mask;
  geno64 <<= 3;
  geno64 |= tmp_g64[1] & trans_mask;
  geno64 <<= 3;
  geno64 |= tmp_g64[0] & trans_mask;
  uint64_t tmp_geno, tmp_mask;
  tmp_geno = geno64 ^ mask_mask;
  geno64 = ~(tmp_geno ^ (tmp_geno << 1));
  tmp_mask = (tmp_geno | (tmp_geno >> 1)) & mask_mask;
  tmp_mask *= 7ull;
  geno64 &= tmp_mask;
  return geno64;
}

void TransposeGeno(const uint8_t *geno8, uint32_t bytes_snp,
                   uint32_t full_bytes, uint32_t indivs_left,
                   uint64_t *geno64) {
  for (uint32_t i = 0; i < full_bytes; ++i) {
    auto g64 = geno64 + i * snps_byte * words_block;
    uint64_t tmp_g64[5];
    ConvertGeno(geno8 + i, bytes_snp, tmp_g64);
    g64[0] = ConvertG64(tmp_g64);
    for (uint32_t j = 1; j < snps_byte; ++j) {
      for (uint32_t k = 0; k < 5; ++k) {
        tmp_g64[k] >>= 2;
      }
      g64[j * words_block] = ConvertG64(tmp_g64);
    }
  }
  if (indivs_left) {
    auto g64 = geno64 + full_bytes * snps_byte * words_block;
    uint64_t tmp_g64[5];
    ConvertGeno(geno8 + full_bytes, bytes_snp, tmp_g64);
    g64[0] = ConvertG64(tmp_g64);
    for (uint32_t j = 1; j < indivs_left; ++j) {
      for (uint32_t k = 0; k < 5; ++k) {
        tmp_g64[k] >>= 2;
      }
      g64[j * words_block] = ConvertG64(tmp_g64);
    }
  }
}

void TransposeLastGeno(const uint8_t *geno8, uint32_t bytes_snp,
                       uint32_t full_bytes, uint32_t indivs_left,
                       uint32_t snps_left, uint64_t *geno64) {
  for (uint32_t i = 0; i < full_bytes; ++i) {
    auto g64 = geno64 + i * snps_byte * words_block;
    uint64_t tmp_g64[5];
    ConvertLastGeno(geno8 + i, bytes_snp, snps_left, tmp_g64);
    g64[0] = ConvertG64(tmp_g64);
    for (uint32_t j = 1; j < snps_byte; ++j) {
      for (uint32_t k = 0; k < 5; ++k) {
        tmp_g64[k] >>= 2;
      }
      g64[j * words_block] = ConvertG64(tmp_g64);
    }
  }
  if (indivs_left) {
    auto g64 = geno64 + full_bytes * snps_byte * words_block;
    uint64_t tmp_g64[4];
    ConvertLastGeno(geno8 + full_bytes, bytes_snp, snps_left, tmp_g64);
    g64[0] = ConvertG64(tmp_g64);
    for (uint32_t j = 1; j < indivs_left; ++j) {
      for (uint32_t k = 0; k < 4; ++k) {
        tmp_g64[k] >>= 2;
      }
      g64[j * words_block] = ConvertG64(tmp_g64);
    }
  }
}

uint64_t CalculateHamming(const uint64_t *geno_1, const uint64_t *geno_2) {
  uint64_t count = 0;
  for (uint32_t i = 0; i < dwords_block; ++i) {
    __m128i g128 = _mm_and_si128(
        _mm_load_si128(reinterpret_cast<const __m128i *>(geno_1 + 2 * i)),
        _mm_load_si128(reinterpret_cast<const __m128i *>(geno_2 + 2 * i)));
    count += _mm_popcnt_u64(_mm_cvtsi128_si64(g128));
    count += _mm_popcnt_u64(_mm_cvtsi128_si64(_mm_srli_si128(g128, 8)));
  }
  return count;
}

void CalcHammingThread(const uint64_t *geno, uint32_t num_indivs,
                       uint32_t num_blocks, uint64_t *matrix, uint32_t begin,
                       uint32_t end) {
  Barrier *br = Barrier::get_instance();
  for (uint32_t l = 0; l < num_blocks; ++l) {
    br->WaitMain();
    for (uint32_t i = begin; i < end; ++i) {
      auto g1 = geno + i * words_block;
      for (uint32_t j = 0; j <= i; ++j) {
        auto g2 = geno + j * words_block;
        matrix[i * num_indivs + j] += CalculateHamming(g1, g2);
      }
    }
    br->WaitThreads();
  }
}
} // namespace
void CalcHammingMatrix(const uint8_t *geno, uint32_t num_indivs,
                       uint32_t num_snps, double *matrix,
                       uint32_t num_threads) {
  uint32_t full_bytes = num_indivs / snps_byte;
  uint32_t indivs_left = num_indivs % snps_byte;
  uint32_t bytes_snp = full_bytes + (indivs_left ? 1 : 0);
  uint32_t full_blocks = num_snps / snps_block;
  uint32_t snps_left = num_snps % snps_block;
  uint32_t num_blocks = full_blocks + (snps_left ? 1 : 0);
  uint64_t *geno64 = new uint64_t[num_indivs * words_block];
  uint64_t *matrix_i = new uint64_t[num_indivs * num_indivs];
  memset(matrix_i, 0, sizeof(uint64_t) * num_indivs * num_indivs);
  Barrier *br = Barrier::get_instance();
  br->SetThreads(num_threads);
  std::vector<std::thread> workers(num_threads);
  double num_jobs = static_cast<double>(num_indivs);
  num_jobs = num_jobs * (num_jobs + 1) / num_threads;
  uint32_t begin = 0;
  for (uint32_t i = 0; i < num_threads - 1; ++i) {
    uint32_t end = static_cast<uint32_t>(
        sqrt(num_jobs + begin * begin + 2 * begin + 1) - 1);
    workers[i] = std::thread(CalcHammingThread, geno64, num_indivs, num_blocks,
                             matrix_i, begin, end);
    begin = end;
  }
  workers[num_threads - 1] =
      std::thread(CalcHammingThread, geno64, num_indivs, num_blocks, matrix_i,
                  begin, num_indivs);
  for (uint32_t i = 0; i < full_blocks; ++i) {
    br->LockMain();
    auto g8 = geno + i * snps_block * bytes_snp;
    for (uint32_t j = 0; j < words_block; ++j) {
      auto g64 = geno64 + j;
      TransposeGeno(g8, bytes_snp, full_bytes, indivs_left, g64);
      g8 += snps_word * bytes_snp;
    }
    br->ReleaseMain();
    br->WaitThreads();
  }
  if (snps_left) {
    br->LockMain();
    uint32_t words_left = snps_left / snps_word;
    auto g8 = geno + full_blocks * snps_block * bytes_snp;
    memset(geno64, 0, sizeof(uint64_t) * words_block * num_indivs);
    for (uint32_t j = 0; j < words_left; ++j) {
      auto g64 = geno64 + j;
      TransposeGeno(g8, bytes_snp, full_bytes, indivs_left, g64);
      g8 += snps_word * bytes_snp;
    }
    snps_left = snps_left % snps_word;
    if (snps_left) {
      auto g64 = geno64 + words_left;
      TransposeLastGeno(g8, bytes_snp, full_bytes, indivs_left, snps_left, g64);
    }
    br->ReleaseMain();
    br->WaitThreads();
  }
  for (auto &&iter : workers) {
    iter.join();
  }
  for (uint32_t i = 0; i < num_indivs; ++i) {
    for (uint32_t j = 0; j < i; ++j) {
      double ibs = static_cast<double>(matrix_i[i * num_indivs + j]);
      ibs /= num_snps;
      matrix[i * num_indivs + j] = ibs;
      matrix[j * num_indivs + i] = ibs;
    }
    double ibs = static_cast<double>(matrix_i[i * num_indivs + i]);
    ibs /= num_snps;
    matrix[i * num_indivs + i] = ibs;
  }
  delete[] geno64;
  delete[] matrix_i;
}
