
#include "poly.h"

#include <immintrin.h>
#include <stdlib.h>

#include "mod.h"
#include "random.h"
#include "utils.h"
#include "xof.h"

void poly_new(poly* m, uint32_t n) { *m = malloc(n * sizeof(int64_t)); }

void poly_free(poly* m) { free(*m); }

// AVX2 doesn't support int64_t <-> double conversion
void vector_to_poly(poly m, const vec* v, uint64_t scale, uint32_t n) {
  // 1/n from IFFT
  double inv_n = (double)scale / n;
  for (uint32_t i = 0; i < n; i++) {
    int64_t mask = ((int64_t)v->rel[i]) >> 63;
    mask = mask & 0x1;
    double round = 0.5 - (double)mask;
    m[i] = (int64_t)(v->rel[i] * inv_n + round);
  }
}

void poly_to_vector(vec* v, const poly m, uint64_t scale, uint32_t n) {
  // To increase the precision, we first compute (1 << int_part)/scale
  // Instead of /scale and *(1 << int_part)
  double scale2 = 1.0 / (double)scale;
  for (uint32_t i = 0; i < n; i++) {
    v->rel[i] = (double)m[i] * scale2;
    v->img[i] = 0;
  }
}

static const uint32_t perm_table[16][8] = {
    // mask     destination 32-bit lanes
    /* 0x0 */ {0, 1, 0, 1, 0, 1, 0, 1},
    /* 0x1 */ {0, 1, 0, 1, 0, 1, 0, 1},
    /* 0x2 */ {2, 3, 0, 1, 0, 1, 0, 1},
    /* 0x3 */ {0, 1, 2, 3, 0, 1, 0, 1},
    /* 0x4 */ {4, 5, 0, 1, 0, 1, 0, 1},
    /* 0x5 */ {0, 1, 4, 5, 0, 1, 0, 1},
    /* 0x6 */ {2, 3, 4, 5, 0, 1, 0, 1},
    /* 0x7 */ {0, 1, 2, 3, 4, 5, 0, 1},
    /* 0x8 */ {6, 7, 0, 1, 0, 1, 0, 1},
    /* 0x9 */ {0, 1, 6, 7, 0, 1, 0, 1},
    /* 0xA */ {2, 3, 6, 7, 0, 1, 0, 1},
    /* 0xB */ {0, 1, 2, 3, 6, 7, 0, 1},
    /* 0xC */ {4, 5, 6, 7, 0, 1, 0, 1},
    /* 0xD */ {0, 1, 4, 5, 6, 7, 0, 1},
    /* 0xE */ {2, 3, 4, 5, 6, 7, 0, 1},
    /* 0xF */ {0, 1, 2, 3, 4, 5, 6, 7}};

void poly_ternary_sampler(poly s, uint32_t n) {
  // Initialize PRNG
  uint8_t seed[SEED_BYTES];
  randombytes(seed, SEED_BYTES);
  xof_state state;
  xof_absorb(&state, seed, SEED_BYTES);
  // Sample
  size_t blocks = (((4 * (n / 4) + 2) / 3) + XOF_RATE - 1) / XOF_RATE;
  size_t buflen = blocks * XOF_RATE;
  uint8_t* buf = malloc(buflen);  // Initial guess
  xof_squeeze(buf, buflen, &state);
  uint32_t i = 0, j = 0;
  __m256i v0123 = _mm256_set_epi64x(3, 2, 1, 0);
  __m256i v0246 = _mm256_set_epi64x(6, 4, 2, 0);
  __m256i v1 = _mm256_set1_epi64x(1);
  __m256i v3 = _mm256_set1_epi64x(3);
  __m256i vn = _mm256_set1_epi64x(n);
  while (i < n) {
    // We're probably almost at the end so just fetch one more block
    if (j == buflen) {
      buflen = XOF_RATE;
      xof_squeeze(buf, buflen, &state);
      j = 0;
    }
    // uint8_t x = buf[j];
    // for (int k = 0; k < 8 / 2 && i < n; k++) {
    //   uint8_t y = (x >> (2 * k)) & 0x3;
    //   if (y < 3) {
    //     int8_t z = (int8_t)(y & 0x1) - (int8_t)(y >> 1);
    //     s[i] = z;
    //     i++;
    //   }
    // }
    __m256i vi = _mm256_set1_epi64x(i);
    vi = _mm256_add_epi64(vi, v0123);
    __m256i maski = _mm256_cmpgt_epi64(vn, vi);
    __m256i vy = _mm256_set1_epi64x(buf[j]);
    vy = _mm256_srlv_epi64(vy, v0246);
    vy = _mm256_and_si256(vy, v3);
    __m256i masky = _mm256_cmpgt_epi64(v3, vy);
    __m256i mask = _mm256_and_si256(maski, masky);
    __m256i vz = _mm256_and_si256(vy, v1);
    __m256i vz2 = _mm256_srlv_epi64(vy, v1);
    vz = _mm256_sub_epi64(vz, vz2);
    uint32_t m = (uint32_t)_mm256_movemask_pd(_mm256_castsi256_pd(mask));
    uint32_t popcnt = (uint32_t)__builtin_popcount(m);
    __m256i perm = _mm256_loadu_si256((__m256i*)perm_table[m]);
    vz = _mm256_permutevar8x32_epi32(vz, perm);
    __m256i store_mask =
        _mm256_set_epi64x(-(int64_t)(popcnt >= 4), -(int64_t)(popcnt >= 3),
                          -(int64_t)(popcnt >= 2), -(int64_t)(popcnt >= 1));
    _mm256_maskstore_epi64((long long*)&s[i], store_mask, vz);
    i += popcnt;
    j++;
  }
  free(buf);
}

// AVX2 doesn't support int64_t multiplication
void poly_uniform_sampler(poly a, uint64_t q, int64_t qRq, uint32_t n) {
  uint8_t seed[SEED_BYTES];
  randombytes(seed, SEED_BYTES);
  // Initialize PRNG
  xof_state state;
  xof_absorb(&state, seed, SEED_BYTES);
  // Sample
  uint64_t bound = UINT64_C(-1) / q;
  bound *= q;
  // Normally q < 2^50, then the rejection probability (2^64-1-kq)/(2^64-1) is
  // extremely small, so I just fetch n int64_t plus the XOF rate rounding
  size_t blocks = (n * sizeof(int64_t));
  blocks = (blocks + XOF_RATE - 1) / XOF_RATE;
  size_t bufbytes = blocks * XOF_RATE;
  // The rate is multiple of 8 so it's fine
  uint64_t* buf = (uint64_t*)malloc(bufbytes);  // Initial guess
  xof_squeeze((uint8_t*)buf, bufbytes, &state);
  uint32_t i = 0, j = 0;
  while (i < n) {
    // We're probably almost at the end so just fetch one more block
    if (j == bufbytes / sizeof(int64_t)) {
      bufbytes = XOF_RATE;
      xof_squeeze((uint8_t*)buf, bufbytes, &state);
      j = 0;
    }
    uint64_t x = buf[j];
    if (x < bound) {
      int64_t y = mod_ct_exact((int64_t)x, (int64_t)q, qRq);
      a[i] = y;
      i++;
    }
    j++;
  }
  free(buf);
}

// AVX2 doesn't support popcnt
void poly_cbd24_sampler_add(poly e, const poly m, uint32_t n) {
  uint8_t* buf = malloc(6 * n * sizeof(uint8_t));
  uint8_t seed[SEED_BYTES];
  randombytes(seed, SEED_BYTES);
  xof(buf, 6 * n, seed, SEED_BYTES);
  for (uint32_t i = 0; i < n; i++) {
    int b1 = __builtin_popcount(buf[6 * i]);
    b1 += __builtin_popcount(buf[6 * i + 1]);
    b1 += __builtin_popcount(buf[6 * i + 2]);
    int b2 = __builtin_popcount(buf[6 * i + 3]);
    b2 += __builtin_popcount(buf[6 * i + 4]);
    b2 += __builtin_popcount(buf[6 * i + 5]);
    e[i] = b1 - b2 + m[i];
  }
  free(buf);
}

void poly_add(poly a, const poly b, uint64_t q, int64_t qRq, int64_t qinv,
              int64_t qR2n, uint32_t n) {
  // R^2/n from INTT
  for (uint32_t i = 0; i < n; i++) {
    a[i] = mod_mul_ct(a[i], qR2n, (int64_t)q, qinv);
    a[i] = mod_ct_exact(a[i] + b[i], (int64_t)q, qRq);
  }
}

void poly_rsub(poly a, const poly b, uint64_t q, int64_t qRq, int64_t qinv,
               int64_t qR2n, uint32_t n) {
  // R^2/n from INTT
  for (uint32_t i = 0; i < n; i++) {
    a[i] = mod_mul_ct(a[i], qR2n, (int64_t)q, qinv);
    a[i] = mod_ct_exact(b[i] - a[i], (int64_t)q, qRq);
  }
}

void poly_pointwise_mul(poly c, const poly a, const poly b, uint64_t q,
                        int64_t qinv, uint32_t n) {
  for (uint32_t i = 0; i < n; i++) {
    c[i] = mod_mul_ct(a[i], b[i], (int64_t)q, qinv);
  }
}

void poly_ntt(poly a, const poly roots, uint64_t q, int64_t qRq, int64_t qinv,
              uint32_t n) {
  int64_t zeta, t;

  for (uint32_t len = n / 2; len > 0; len >>= 1) {
    uint32_t k = n / (2 * len);
    for (uint32_t start = 0; start < n; start += 2 * len) {
      zeta = roots[k];
      k++;
      for (uint32_t j = start; j < start + len; ++j) {
        t = mod_mul_ct(zeta, a[j + len], (int64_t)q, qinv);
        a[j + len] = mod_ct(a[j] - t, (int64_t)q, qRq);
        a[j] = mod_ct(a[j] + t, (int64_t)q, qRq);
      }
    }
  }
}

void poly_intt(poly a, const poly roots, uint64_t q, int64_t qRq, int64_t qinv,
               uint32_t n) {
  int64_t t, zeta;

  for (uint32_t len = 1; len < n; len <<= 1) {
    uint32_t k = n / (2 * len);
    for (uint32_t start = 0; start < n; start += 2 * len) {
      zeta = roots[k];
      k++;
      for (uint32_t j = start; j < start + len; ++j) {
        t = a[j];
        a[j] = mod_ct(t + a[j + len], (int64_t)q, qRq);
        a[j + len] = t - a[j + len];  // lazy reduction with the multiplication
        a[j + len] = mod_mul_ct(zeta, a[j + len], (int64_t)q, qinv);
      }
    }
  }

  // R^2/n multiplication performed in poly_add/poly_rsub
}

void poly_gen_roots(poly roots, uint32_t n, int64_t q, int64_t qR,
                    int64_t qinv) {
  // Find one primitive root (n > 2)
  uint32_t log2_n = floor_log2(n);
  int64_t root0 = find_smallest_2nth_root(n, q);
  root0 = mod_mul(root0, qR, q);  // montgomery factor
  roots[0] = 1;
  roots[n / 2] = root0;
  for (uint32_t i = 2; i < n; i++) {
    roots[bit_reverse32(i, log2_n)] =
        mod_mul_ct_exact(roots[bit_reverse32(i - 1, log2_n)], root0, q, qinv);
  }
}

void poly_gen_inverse_roots(poly roots, uint32_t n, int64_t q, int64_t qinv) {
  // Find one primitive root (n > 2)
  uint32_t log2_n = floor_log2(n);
  int64_t etan = -1;             // x^n
  int64_t root0 = roots[n / 2];  // Previously poly_gen_roots output
  roots[0] = 1;
  etan = mod_mul(etan, root0, q);
  roots[n - 1] = etan;  // bit_reverse32(n - 1, log2_n)
  for (uint32_t i = 2; i < n; i++) {
    etan = mod_mul_ct_exact(etan, root0, q, qinv);
    roots[bit_reverse32(n - i, log2_n)] = etan;
  }
}
