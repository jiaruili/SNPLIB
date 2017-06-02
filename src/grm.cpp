#include "grm.h"

namespace {
const uint32_t snps_byte = 4;
const uint32_t snps_word = 20;
const uint64_t trans_mask = 0x0003000300030003ull; // 0003 = 0000000000000011
const uint64_t mask_mask = 0x1249124912491249ull;  // 1249 = 0001001001001001

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
  return geno64;
}

void TransposeGeno(const uint8_t *geno8, uint32_t bytes_snp,
                   uint32_t full_bytes, uint32_t indivs_left,
                   uint64_t *geno64) {
  for (uint32_t i = 0; i < full_bytes; ++i) {
    auto g64 = geno64 + i * snps_byte;
    uint64_t tmp_g64[5];
    ConvertGeno(geno8 + i, bytes_snp, tmp_g64);
    g64[0] = ConvertG64(tmp_g64);
    for (uint32_t j = 1; j < snps_byte; ++j) {
      for (uint32_t k = 0; k < 5; ++k) {
        tmp_g64[k] >>= 2;
      }
      g64[j] = ConvertG64(tmp_g64);
    }
  }
  if (indivs_left) {
    auto g64 = geno64 + full_bytes * snps_byte;
    uint64_t tmp_g64[5];
    ConvertGeno(geno8 + full_bytes, bytes_snp, tmp_g64);
    g64[0] = ConvertG64(tmp_g64);
    for (uint32_t j = 1; j < indivs_left; ++j) {
      for (uint32_t k = 0; k < 5; ++k) {
        tmp_g64[k] >>= 2;
      }
      g64[j] = ConvertG64(tmp_g64);
    }
  }
}

void TransposeLastGeno(const uint8_t *geno8, uint32_t bytes_snp,
                       uint32_t full_bytes, uint32_t indivs_left,
                       uint32_t snps_left, uint64_t *geno64) {
  for (uint32_t i = 0; i < full_bytes; ++i) {
    auto g64 = geno64 + i * snps_byte;
    uint64_t tmp_g64[5];
    ConvertLastGeno(geno8 + i, bytes_snp, snps_left, tmp_g64);
    g64[0] = ConvertG64(tmp_g64);
    for (uint32_t j = 1; j < snps_byte; ++j) {
      for (uint32_t k = 0; k < 5; ++k) {
        tmp_g64[k] >>= 2;
      }
      g64[j] = ConvertG64(tmp_g64);
    }
  }
  if (indivs_left) {
    auto g64 = geno64 + full_bytes * snps_byte;
    uint64_t tmp_g64[4];
    ConvertLastGeno(geno8 + full_bytes, bytes_snp, snps_left, tmp_g64);
    g64[0] = ConvertG64(tmp_g64);
    for (uint32_t j = 1; j < indivs_left; ++j) {
      for (uint32_t k = 0; k < 4; ++k) {
        tmp_g64[k] >>= 2;
      }
      g64[j] = ConvertG64(tmp_g64);
    }
  }
}

void MaskGeno(const uint64_t *geno64, uint32_t num_indivs, uint64_t *mask64) {
  for (uint32_t i = 0; i < num_indivs; ++i) {
    uint64_t tmp_geno = ~(geno64[i] ^ mask_mask);
    tmp_geno = (tmp_geno & (tmp_geno >> 1)) & mask_mask;
    mask64[i] = tmp_geno * 7ull;
  }
}

void CalcSnpCorr(double table[8], double af) {
  table[0] = 2.0 * af / (1.0 - af);       // 00 00
  table[1] = 0.0;                         // 00 01
  table[2] = 1.0 / (1.0 - af) - 2;        // 00 10
  table[3] = -2.0;                        // 00 11
  table[4] = 0.5 / af / (1.0 - af) - 2.0; // 10 10
  table[5] = 1.0 / af - 2.0;              // 10 11
  table[6] = 2.0 / af - 2.0;              // 11 11
  table[7] = 0.0;                         // 01 XX
}

void UpdateLookupTable(const double af[snps_word],
                       double lookup_table[131072]) {
  for (uint32_t i = 0; i < 4; ++i) {
    auto lkt = lookup_table + i * 32768;
    double table[40];
    for (uint32_t j = 0; j < 5; j++) {
      CalcSnpCorr(table + 8 * j, af[5 * i + j]);
    }
    for (uint32_t j = 0; j < 32768; ++j) {
      lkt[j] = table[j & 7u] + table[8 + (j >> 3 & 7u)] +
               table[16 + (j >> 6 & 7u)] + table[24 + (j >> 9 & 7u)] +
               table[32 + (j >> 12 & 7u)];
    }
  }
}

void UpdateLookupTableLast(const double *af, uint32_t snps_left,
                           double lookup_table[131072]) {
  uint32_t shorts_left = snps_left / 5;
  uint32_t snps_remain = snps_left % 5;
  for (uint32_t i = 0; i < shorts_left; ++i) {
    auto lkt = lookup_table + i * 32768;
    double table[40];
    for (uint32_t j = 0; j < 5; j++) {
      CalcSnpCorr(table + 8 * j, af[5 * i + j]);
    }
    for (uint32_t j = 0; j < 32768; ++j) {
      lkt[j] = table[j & 7u] + table[8 + (j >> 3 & 7u)] +
               table[16 + (j >> 6 & 7u)] + table[24 + (j >> 9 & 7u)] +
               table[32 + (j >> 12 & 7u)];
    }
  }
  if (snps_remain) {
    auto lkt = lookup_table + shorts_left * 32768;
    double table[40]{0.0};
    for (uint32_t j = 0; j < snps_remain; j++) {
      CalcSnpCorr(table + 8 * j, af[5 * shorts_left + j]);
    }
    for (uint32_t j = 0; j < 32768; ++j) {
      lkt[j] = table[j & 7u] + table[8 + (j >> 3 & 7u)] +
               table[16 + (j >> 6 & 7u)] + table[24 + (j >> 9 & 7u)] +
               table[32 + (j >> 12 & 7u)];
    }
  }
}

void CalcGrmThread(const uint64_t *geno, const uint64_t *mask,
                   const double lookup_table[131072], uint32_t num_indivs,
                   uint32_t num_words, double *matrix, uint32_t begin,
                   uint32_t end) {
  Barrier *br = Barrier::get_instance();
  for (uint32_t l = 0; l < num_words; ++l) {
    br->WaitMain();
    for (uint32_t i = begin; i < end; ++i) {
      if (!mask[i]) {
        for (uint32_t j = 0; j <= i; ++j) {
          uint64_t g = (geno[i] + geno[j]) | mask[j];
          matrix[i * num_indivs + j] +=
              lookup_table[g & 0x7fffu] +
              lookup_table[(g >> 16 & 0x7fffu) + 0x8000u] +
              lookup_table[(g >> 32 & 0x7fffu) + 0x10000u] +
              lookup_table[(g >> 48 & 0x7fffu) + 0x18000u];
        }
      } else {
        for (uint32_t j = 0; j <= i; ++j) {
          uint64_t g = (geno[i] + geno[j]) | (mask[i] | mask[j]);
          matrix[i * num_indivs + j] +=
              lookup_table[g & 0x7fffu] +
              lookup_table[(g >> 16 & 0x7fffu) + 0x8000u] +
              lookup_table[(g >> 32 & 0x7fffu) + 0x10000u] +
              lookup_table[(g >> 48 & 0x7fffu) + 0x18000u];
        }
      }
    }
    br->WaitThreads();
  }
}
} // namespace

void CalcGrmMatrix(const uint8_t *geno, const double *af, uint32_t num_indivs,
                   uint32_t num_snps, double *matrix, uint32_t num_threads) {
  uint32_t full_bytes = num_indivs / snps_byte;
  uint32_t indivs_left = num_indivs % snps_byte;
  uint32_t bytes_snp = full_bytes + (indivs_left ? 1 : 0);
  uint32_t full_words = num_snps / snps_word;
  uint32_t snps_left = num_snps % snps_word;
  uint32_t num_words = full_words + (snps_left ? 1 : 0);
  uint64_t *geno64 = new uint64_t[num_indivs];
  uint64_t *mask64 = new uint64_t[num_indivs];
  double lookup_table[131072];
  Barrier *br = Barrier::get_instance();
  br->SetThreads(num_threads);
  std::vector<std::thread> workers(num_threads);
  double num_jobs = static_cast<double>(num_indivs);
  num_jobs = num_jobs * (num_jobs + 1) / num_threads;
  uint32_t begin = 0;
  for (uint32_t i = 0; i < num_threads - 1; ++i) {
    uint32_t end = static_cast<uint32_t>(
        sqrt(num_jobs + begin * begin + 2 * begin + 1) - 1);
    workers[i] = std::thread(CalcGrmThread, geno64, mask64, lookup_table,
                             num_indivs, num_words, matrix, begin, end);
    begin = end;
  }
  workers[num_threads - 1] =
      std::thread(CalcGrmThread, geno64, mask64, lookup_table, num_indivs,
                  num_words, matrix, begin, num_indivs);
  for (uint32_t i = 0; i < full_words; ++i) {
    br->LockMain();
    auto g8 = geno + i * snps_word * bytes_snp;
    TransposeGeno(g8, bytes_snp, full_bytes, indivs_left, geno64);
    MaskGeno(geno64, num_indivs, mask64);
    auto tmp_af = af + i * snps_word;
    UpdateLookupTable(tmp_af, lookup_table);
    br->ReleaseMain();
    br->WaitThreads();
  }
  if (snps_left) {
    br->LockMain();
    auto g8 = geno + full_words * snps_word * bytes_snp;
    TransposeLastGeno(g8, bytes_snp, full_bytes, indivs_left, snps_left,
                      geno64);
    MaskGeno(geno64, num_indivs, mask64);
    auto tmp_af = af + full_words * snps_word;
    UpdateLookupTableLast(tmp_af, snps_left, lookup_table);
    br->ReleaseMain();
    br->WaitThreads();
  }
  for (auto &&iter : workers) {
    iter.join();
  }
  for (uint32_t i = 0; i < num_indivs; ++i) {
    for (uint32_t j = 0; j < i; ++j) {
      matrix[i * num_indivs + j] /= num_snps;
      matrix[j * num_indivs + i] = matrix[i * num_indivs + j];
    }
    matrix[i * num_indivs + i] /= num_snps;
  }
  delete[] geno64;
  delete[] mask64;
}
