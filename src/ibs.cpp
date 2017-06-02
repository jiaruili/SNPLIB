#include "ibs.h"

namespace {
const uint32_t snps_byte = 4;
const uint32_t snps_word = 32;
const uint32_t words_block = 32;
const uint64_t trans_mask = 0x0303030303030303ull; // 0303 = 0000001100000011
const uint64_t mask_mask = 0x5555555555555555ull;  // 5555 = 0101010101010101
const auto snps_block = snps_word * words_block;
const auto dwords_block = words_block / 2;

void ConvertGeno(const uint8_t *geno, uint32_t bytes_snp, uint64_t tmp_g64[4]) {
  uint64_t tmp_g8[snps_word];
  for (uint32_t i = 0; i < snps_word; ++i) {
    tmp_g8[i] = static_cast<uint64_t>(geno[i * bytes_snp]);
  }
  for (uint32_t i = 0; i < 4; ++i) {
    tmp_g64[i] = tmp_g8[i] + (tmp_g8[4 + i] << 8) + (tmp_g8[8 + i] << 16) +
                 (tmp_g8[12 + i] << 24) + (tmp_g8[16 + i] << 32) +
                 (tmp_g8[20 + i] << 40) + (tmp_g8[24 + i] << 48) +
                 (tmp_g8[28 + i] << 56);
  }
}

void ConvertLastGeno(const uint8_t *geno, uint32_t bytes_snp,
                     uint32_t snps_left, uint64_t tmp_g64[4]) {
  uint64_t tmp_g8[snps_word];
  for (uint32_t i = 0; i < snps_left; ++i) {
    tmp_g8[i] = static_cast<uint64_t>(geno[i * bytes_snp]);
  }
  for (uint32_t i = snps_left; i < snps_word; ++i) {
    tmp_g8[i] = 0x55ull;
  }
  for (uint32_t i = 0; i < 4; ++i) {
    tmp_g64[i] = tmp_g8[i] + (tmp_g8[4 + i] << 8) + (tmp_g8[8 + i] << 16) +
                 (tmp_g8[12 + i] << 24) + (tmp_g8[16 + i] << 32) +
                 (tmp_g8[20 + i] << 40) + (tmp_g8[24 + i] << 48) +
                 (tmp_g8[28 + i] << 56);
  }
}

uint64_t ConvertG64(const uint64_t tmp_g64[4]) {
  uint64_t geno64;
  geno64 = tmp_g64[3] & trans_mask;
  geno64 <<= 2;
  geno64 |= tmp_g64[2] & trans_mask;
  geno64 <<= 2;
  geno64 |= tmp_g64[1] & trans_mask;
  geno64 <<= 2;
  geno64 |= tmp_g64[0] & trans_mask;
  return geno64;
}

void TransposeGeno(const uint8_t *geno8, uint32_t bytes_snp,
                   uint32_t full_bytes, uint32_t indivs_left,
                   uint64_t *geno64) {
  for (uint32_t i = 0; i < full_bytes; ++i) {
    auto g64 = geno64 + i * snps_byte * words_block;
    uint64_t tmp_g64[4];
    ConvertGeno(geno8 + i, bytes_snp, tmp_g64);
    g64[0] = ConvertG64(tmp_g64);
    for (uint32_t j = 1; j < snps_byte; ++j) {
      for (uint32_t k = 0; k < 4; ++k) {
        tmp_g64[k] >>= 2;
      }
      g64[j * words_block] = ConvertG64(tmp_g64);
    }
  }
  if (indivs_left) {
    auto g64 = geno64 + full_bytes * snps_byte * words_block;
    uint64_t tmp_g64[4];
    ConvertGeno(geno8 + full_bytes, bytes_snp, tmp_g64);
    g64[0] = ConvertG64(tmp_g64);
    for (uint32_t j = 1; j < indivs_left; ++j) {
      for (uint32_t k = 0; k < 4; ++k) {
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
    uint64_t tmp_g64[4];
    ConvertLastGeno(geno8 + i, bytes_snp, snps_left, tmp_g64);
    g64[0] = ConvertG64(tmp_g64);
    for (uint32_t j = 1; j < snps_byte; ++j) {
      for (uint32_t k = 0; k < 4; ++k) {
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

void MaskGeno(const uint64_t *geno64, uint32_t num_indivs, uint64_t *mask64) {
  for (uint32_t i = 0; i < num_indivs * words_block; ++i) {
    uint64_t tmp_geno = geno64[i] ^ mask_mask;
    tmp_geno = (tmp_geno | (tmp_geno >> 1)) & mask_mask;
    mask64[i] = tmp_geno * 3ull;
  }
}

uint64_t CalcIbs(const uint64_t *geno_1, const uint64_t *geno_2,
                 const uint64_t *mask_1, const uint64_t *mask_2) {
  uint64_t count = 0;
  for (uint32_t l = 0; l < dwords_block; ++l) {
    __m128i g = _mm_xor_si128(
        _mm_load_si128(reinterpret_cast<const __m128i *>(geno_1 + 2 * l)),
        _mm_load_si128(reinterpret_cast<const __m128i *>(geno_2 + 2 * l)));
    __m128i m = _mm_and_si128(
        _mm_load_si128(reinterpret_cast<const __m128i *>(mask_1 + 2 * l)),
        _mm_load_si128(reinterpret_cast<const __m128i *>(mask_2 + 2 * l)));
    g = _mm_andnot_si128(g, m);
    count += _mm_popcnt_u64(static_cast<uint64_t>(_mm_cvtsi128_si64(g)));
    count += _mm_popcnt_u64(
        static_cast<uint64_t>(_mm_cvtsi128_si64(_mm_srli_si128(g, 8))));
  }
  return count;
}

void CalcIbsThread(const uint64_t *geno, const uint64_t *mask,
                   uint32_t num_indivs, uint32_t num_blocks, uint64_t *matrix,
                   uint32_t begin, uint32_t end) {
  Barrier *br = Barrier::get_instance();
  for (uint32_t l = 0; l < num_blocks; ++l) {
    br->WaitMain();
    for (uint32_t i = begin; i < end; ++i) {
      auto g1 = geno + i * words_block;
      auto m1 = mask + i * words_block;
      for (uint32_t j = 0; j <= i; ++j) {
        auto g2 = geno + j * words_block;
        auto m2 = mask + j * words_block;
        matrix[i * num_indivs + j] += CalcIbs(g1, g2, m1, m2);
      }
    }
    br->WaitThreads();
  }
}

void CalcConnectThread(const uint64_t *geno_src, const uint64_t *mask_src,
                       const uint64_t *geno_dest, const uint64_t *mask_dest,
                       uint32_t num_indivs_src, uint32_t num_indivs_dest,
                       uint32_t num_blocks, uint64_t *connect) {
  Barrier *br = Barrier::get_instance();
  for (uint32_t l = 0; l < num_blocks; ++l) {
    br->WaitMain();
    for (uint32_t i = 0; i < num_indivs_dest; ++i) {
      auto g1 = geno_dest + i * words_block;
      auto m1 = mask_dest + i * words_block;
      for (uint32_t j = 0; j < num_indivs_src; j++) {
        auto g2 = geno_src + j * words_block;
        auto m2 = mask_src + j * words_block;
        connect[i] += CalcIbs(g1, g2, m1, m2);
      }
    }
    br->WaitThreads();
  }
}
} // namespace
void CalcIbsMatrix(const uint8_t *geno, uint32_t num_indivs, uint32_t num_snps,
                   double *matrix, uint32_t num_threads) {
  uint32_t full_bytes = num_indivs / snps_byte;
  uint32_t indivs_left = num_indivs % snps_byte;
  uint32_t bytes_snp = full_bytes + (indivs_left ? 1 : 0);
  uint32_t full_blocks = num_snps / snps_block;
  uint32_t snps_left = num_snps % snps_block;
  uint32_t num_blocks = full_blocks + (snps_left ? 1 : 0);
  uint64_t *geno64 = new uint64_t[num_indivs * words_block];
  uint64_t *mask64 = new uint64_t[num_indivs * words_block];
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
    workers[i] = std::thread(CalcIbsThread, geno64, mask64, num_indivs,
                             num_blocks, matrix_i, begin, end);
    begin = end;
  }
  workers[num_threads - 1] =
      std::thread(CalcIbsThread, geno64, mask64, num_indivs, num_blocks,
                  matrix_i, begin, num_indivs);
  for (uint32_t i = 0; i < full_blocks; ++i) {
    br->LockMain();
    auto g8 = geno + i * snps_block * bytes_snp;
    for (uint32_t j = 0; j < words_block; ++j) {
      auto g64 = geno64 + j;
      TransposeGeno(g8, bytes_snp, full_bytes, indivs_left, g64);
      g8 += snps_word * bytes_snp;
    }
    MaskGeno(geno64, num_indivs, mask64);
    br->ReleaseMain();
    br->WaitThreads();
  }
  if (snps_left) {
    br->LockMain();
    uint32_t words_left = snps_left / snps_word;
    auto g8 = geno + full_blocks * snps_block * bytes_snp;
    memset(geno64, 0x55u, sizeof(uint64_t) * words_block * num_indivs);
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
    MaskGeno(geno64, num_indivs, mask64);
    br->ReleaseMain();
    br->WaitThreads();
  }
  for (auto &&iter : workers) {
    iter.join();
  }
  for (uint32_t i = 0; i < num_indivs; ++i) {
    for (uint32_t j = 0; j < i; ++j) {
      double ibs = static_cast<double>(matrix_i[i * num_indivs + j]);
      ibs /= 2 * num_snps;
      matrix[i * num_indivs + j] = ibs;
      matrix[j * num_indivs + i] = ibs;
    }
    double ibs = static_cast<double>(matrix_i[i * num_indivs + i]);
    ibs /= 2 * num_snps;
    matrix[i * num_indivs + i] = ibs;
  }
  delete[] geno64;
  delete[] mask64;
  delete[] matrix_i;
}

void CalcIbsConnection(const uint8_t *geno_src, uint32_t num_indivs_src,
                       const uint8_t *geno_dest, uint32_t num_indivs_dest,
                       uint32_t num_snps, double *connection,
                       uint32_t num_threads) {
  uint32_t full_bytes_src = num_indivs_src / snps_byte;
  uint32_t indivs_left_src = num_indivs_src % snps_byte;
  uint32_t bytes_snp_src = full_bytes_src + (indivs_left_src ? 1 : 0);
  uint32_t full_bytes_dest = num_indivs_dest / snps_byte;
  uint32_t indivs_left_dest = num_indivs_dest % snps_byte;
  uint32_t bytes_snp_dest = full_bytes_dest + (indivs_left_dest ? 1 : 0);
  uint32_t full_blocks = num_snps / snps_block;
  uint32_t snps_left = num_snps % snps_block;
  uint32_t num_blocks = full_blocks + (snps_left ? 1 : 0);
  uint64_t *geno64_src = new uint64_t[num_indivs_src * words_block];
  uint64_t *mask64_src = new uint64_t[num_indivs_src * words_block];
  uint64_t *geno64_dest = new uint64_t[num_indivs_dest * words_block];
  uint64_t *mask64_dest = new uint64_t[num_indivs_dest * words_block];
  uint64_t *connect = new uint64_t[num_indivs_dest];
  memset(connect, 0, sizeof(uint64_t) * num_indivs_dest);
  Barrier *br = Barrier::get_instance();
  br->SetThreads(num_threads);
  std::vector<std::thread> workers(num_threads);
  uint32_t indivs_job = num_indivs_dest / num_threads + 1;
  uint32_t indivs_left = num_indivs_dest % num_threads;
  uint64_t *g64 = geno64_dest;
  uint64_t *m64 = mask64_dest;
  for (uint32_t i = 0; i < indivs_left; ++i) {
    workers[i] =
        std::thread(CalcConnectThread, geno64_src, mask64_src, g64, m64,
                    num_indivs_src, indivs_job, num_blocks, connect);
    g64 += indivs_job * words_block;
    m64 += indivs_job * words_block;
    connect += indivs_job;
  }
  --indivs_job;
  for (uint32_t i = indivs_left; i < num_threads; ++i) {
    workers[i] =
        std::thread(CalcConnectThread, geno64_src, mask64_src, g64, m64,
                    num_indivs_src, indivs_job, num_blocks, connect);
    g64 += indivs_job * words_block;
    m64 += indivs_job * words_block;
    connect += indivs_job;
  }
  for (uint32_t i = 0; i < full_blocks; ++i) {
    br->LockMain();
    auto g8 = geno_src + i * snps_block * bytes_snp_src;
    for (uint32_t j = 0; j < words_block; ++j) {
      auto g64 = geno64_src + j;
      TransposeGeno(g8, bytes_snp_src, full_bytes_src, indivs_left_src, g64);
      g8 += snps_word * bytes_snp_src;
    }
    g8 = geno_dest + i * snps_block * bytes_snp_dest;
    for (uint32_t j = 0; j < words_block; ++j) {
      auto g64 = geno64_dest + j;
      TransposeGeno(g8, bytes_snp_dest, full_bytes_dest, indivs_left_dest, g64);
      g8 += snps_word * bytes_snp_dest;
    }
    MaskGeno(geno64_src, num_indivs_src, mask64_src);
    MaskGeno(geno64_dest, num_indivs_dest, mask64_dest);
    br->ReleaseMain();
    br->WaitThreads();
  }
  if (snps_left) {
    br->LockMain();
    uint32_t words_left = snps_left / snps_word;
    auto g8 = geno_src + full_blocks * snps_block * bytes_snp_src;
    memset(geno64_src, 0x55u, sizeof(uint64_t) * words_block * num_indivs_src);
    for (uint32_t j = 0; j < words_left; ++j) {
      auto g64 = geno64_src + j;
      TransposeGeno(g8, bytes_snp_src, full_bytes_src, indivs_left_src, g64);
      g8 += snps_word * bytes_snp_src;
    }
    g8 = geno_dest + full_blocks * snps_block * bytes_snp_dest;
    memset(geno64_dest, 0x55u,
           sizeof(uint64_t) * words_block * num_indivs_dest);
    for (uint32_t j = 0; j < words_left; ++j) {
      auto g64 = geno64_dest + j;
      TransposeGeno(g8, bytes_snp_dest, full_bytes_dest, indivs_left_dest, g64);
      g8 += snps_word * bytes_snp_dest;
    }
    snps_left = snps_left % snps_word;
    if (snps_left) {
      auto g64 = geno64_src + words_left;
      TransposeLastGeno(g8, bytes_snp_src, full_bytes_src, indivs_left_src,
                        snps_left, g64);
      g64 = geno64_dest + words_left;
      TransposeLastGeno(g8, bytes_snp_dest, full_bytes_dest, indivs_left_dest,
                        snps_left, g64);
    }
    MaskGeno(geno64_src, num_indivs_src, mask64_src);
    MaskGeno(geno64_dest, num_indivs_dest, mask64_dest);
    br->ReleaseMain();
    br->WaitThreads();
  }
  for (auto &&iter : workers) {
    iter.join();
  }
  for (uint32_t i = 0; i < num_indivs_dest; ++i) {
    connection[i] = static_cast<double>(connect[i]) / 2 / num_snps;
  }
  delete[] geno64_src;
  delete[] mask64_src;
  delete[] geno64_dest;
  delete[] mask64_dest;
}
