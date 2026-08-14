
#include "poly.h"

#include <immintrin.h>
#include <stdlib.h>

#include "mod.h"
#include "random.h"
#include "utils.h"
#include "xof.h"

void poly_new(poly_t* m, uint32_t n) { *m = malloc(n * sizeof(int64_t)); }

void poly_free(poly_t* m) { free(*m); }

// AVX2 doesn't support int64_t <-> double conversion
void vector_to_poly(poly_t m, const vec_t* v, uint64_t scale, uint32_t n) {
  // 1/n from IFFT
  double inv_n = (double)scale / n;
  for (uint32_t i = 0; i < n; i++) {
    int64_t mask = ((int64_t)v->rel[i]) >> 63;
    mask = mask & 0x1;
    double round = 0.5 - (double)mask;
    m[i] = (int64_t)(v->rel[i] * inv_n + round);
  }
}

void poly_to_vector(vec_t* v, const poly_t m, uint64_t scale, uint32_t n) {
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

void poly_ternary_sampler(poly_t s, uint32_t n) {
  // Initialize PRNG
  uint8_t seed[SEED_BYTES];
  randombytes(seed, SEED_BYTES);
  xof_state_t state;
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
void poly_uniform_sampler(poly_t a, uint64_t q, int64_t qRq, uint32_t n) {
  uint8_t seed[SEED_BYTES];
  randombytes(seed, SEED_BYTES);
  // Initialize PRNG
  xof_state_t state;
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
void poly_cbd24_sampler_add(poly_t e, const poly_t m, uint32_t n) {
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

void poly_add(poly_t a, const poly_t b, uint64_t q, int64_t qRq, int64_t qinv,
              int64_t qR2n, uint32_t n) {
  // R^2/n from INTT
  for (uint32_t i = 0; i < n; i++) {
    a[i] = mod_mul_ct(a[i], qR2n, (int64_t)q, qinv);
    a[i] = mod_ct_exact(a[i] + b[i], (int64_t)q, qRq);
  }
}

void poly_rsub(poly_t a, const poly_t b, uint64_t q, int64_t qRq, int64_t qinv,
               int64_t qR2n, uint32_t n) {
  // R^2/n from INTT
  for (uint32_t i = 0; i < n; i++) {
    a[i] = mod_mul_ct(a[i], qR2n, (int64_t)q, qinv);
    a[i] = mod_ct_exact(b[i] - a[i], (int64_t)q, qRq);
  }
}

void poly_pointwise_mul(poly_t c, const poly_t a, const poly_t b, uint64_t q,
                        int64_t qinv, uint32_t n) {
  for (uint32_t i = 0; i < n; i++) {
    c[i] = mod_mul_ct(a[i], b[i], (int64_t)q, qinv);
  }
}

#define CT_BUTTERFLY(aj, ak, zeta, t)                 \
  do {                                                \
    (t) = mod_mul_ct((zeta), (ak), (int64_t)q, qinv); \
    (ak) = mod_ct((aj) - (t), (int64_t)q, qRq);       \
    (aj) = mod_ct((aj) + (t), (int64_t)q, qRq);       \
  } while (0)

void poly_ntt(poly_t a, const poly_t roots, uint64_t q, int64_t qRq,
              int64_t qinv, uint32_t n) {
  int64_t zeta, t;

  uint32_t log2_n = floor_log2(n / 16);
  uint32_t n3_rem = log2_n % 3;
  if (n3_rem == 1) {  // Just First layer
    // Roots
    int64_t w = roots[1];
    for (uint32_t j = 0; j < n / 2; j++) {
      // Load
      int64_t aj0 = a[j + 0 * (n / 2)];
      int64_t aj1 = a[j + 1 * (n / 2)];
      // First Layer
      CT_BUTTERFLY(aj0, aj1, w, t);
      // Store
      a[j + 0 * (n / 2)] = aj0;
      a[j + 1 * (n / 2)] = aj1;
    }
  } else if (n3_rem == 2) {  // First and Second layer
    // Roots
    int64_t w1 = roots[1];
    int64_t w20 = roots[2];
    int64_t w21 = roots[3];
    for (uint32_t j = 0; j < n / 4; j++) {
      // Load
      int64_t aj0 = a[j + 0 * (n / 4)];
      int64_t aj1 = a[j + 1 * (n / 4)];
      int64_t aj2 = a[j + 2 * (n / 4)];
      int64_t aj3 = a[j + 3 * (n / 4)];
      // First layer (len)
      CT_BUTTERFLY(aj0, aj2, w1, t);
      CT_BUTTERFLY(aj1, aj3, w1, t);
      // Second layer (len/2)
      CT_BUTTERFLY(aj0, aj1, w20, t);
      CT_BUTTERFLY(aj2, aj3, w21, t);
      // Store
      a[j + 0 * (n / 4)] = aj0;
      a[j + 1 * (n / 4)] = aj1;
      a[j + 2 * (n / 4)] = aj2;
      a[j + 3 * (n / 4)] = aj3;
    }
  }

  uint32_t n3 = 1 << (log2_n + 4 - n3_rem);
  // Merge 3 layers using 16+ registers
  for (uint32_t len = n3 / 2; len > 8; len >>= 3) {
    uint32_t k1 = n / (2 * len);
    uint32_t k2 = k1 << 1;
    uint32_t k3 = k1 << 2;
    for (uint32_t start = 0; start < n; start += 2 * len) {
      // Roots
      int64_t w1 = roots[k1];
      k1 += 1;
      int64_t w20 = roots[k2];
      int64_t w21 = roots[k2 + 1];
      k2 += 2;
      int64_t w30 = roots[k3];
      int64_t w31 = roots[k3 + 1];
      int64_t w32 = roots[k3 + 2];
      int64_t w33 = roots[k3 + 3];
      k3 += 4;
      for (uint32_t j = start; j < start + len / 4; j++) {
        // Load
        int64_t aj0 = a[j + 0 * (len / 4)];
        int64_t aj1 = a[j + 1 * (len / 4)];
        int64_t aj2 = a[j + 2 * (len / 4)];
        int64_t aj3 = a[j + 3 * (len / 4)];
        int64_t aj4 = a[j + 4 * (len / 4)];
        int64_t aj5 = a[j + 5 * (len / 4)];
        int64_t aj6 = a[j + 6 * (len / 4)];
        int64_t aj7 = a[j + 7 * (len / 4)];
        // First layer (len), one root
        CT_BUTTERFLY(aj0, aj4, w1, t);
        CT_BUTTERFLY(aj1, aj5, w1, t);
        CT_BUTTERFLY(aj2, aj6, w1, t);
        CT_BUTTERFLY(aj3, aj7, w1, t);
        // Second layer (len/2), one root per half
        CT_BUTTERFLY(aj0, aj2, w20, t);
        CT_BUTTERFLY(aj1, aj3, w20, t);
        CT_BUTTERFLY(aj4, aj6, w21, t);
        CT_BUTTERFLY(aj5, aj7, w21, t);
        // Third layer (len/4), one root per quarter
        CT_BUTTERFLY(aj0, aj1, w30, t);
        CT_BUTTERFLY(aj2, aj3, w31, t);
        CT_BUTTERFLY(aj4, aj5, w32, t);
        CT_BUTTERFLY(aj6, aj7, w33, t);
        // Store
        a[j + 0 * (len / 4)] = aj0;
        a[j + 1 * (len / 4)] = aj1;
        a[j + 2 * (len / 4)] = aj2;
        a[j + 3 * (len / 4)] = aj3;
        a[j + 4 * (len / 4)] = aj4;
        a[j + 5 * (len / 4)] = aj5;
        a[j + 6 * (len / 4)] = aj6;
        a[j + 7 * (len / 4)] = aj7;
      }
    }
  }

  // Merge last 4 layers using 16+ registers
  uint32_t k1 = n / 16;
  uint32_t k2 = n / 8;
  uint32_t k3 = n / 4;
  uint32_t k4 = n / 2;
  for (uint32_t j = 0; j < n; j += 16) {
    // Load
    int64_t aj0 = a[j + 0];
    int64_t aj1 = a[j + 1];
    int64_t aj2 = a[j + 2];
    int64_t aj3 = a[j + 3];
    int64_t aj4 = a[j + 4];
    int64_t aj5 = a[j + 5];
    int64_t aj6 = a[j + 6];
    int64_t aj7 = a[j + 7];
    int64_t aj8 = a[j + 8];
    int64_t aj9 = a[j + 9];
    int64_t aj10 = a[j + 10];
    int64_t aj11 = a[j + 11];
    int64_t aj12 = a[j + 12];
    int64_t aj13 = a[j + 13];
    int64_t aj14 = a[j + 14];
    int64_t aj15 = a[j + 15];
    // First layer (len = 8), one root
    zeta = roots[k1];
    k1 += 1;
    CT_BUTTERFLY(aj0, aj8, zeta, t);
    CT_BUTTERFLY(aj1, aj9, zeta, t);
    CT_BUTTERFLY(aj2, aj10, zeta, t);
    CT_BUTTERFLY(aj3, aj11, zeta, t);
    CT_BUTTERFLY(aj4, aj12, zeta, t);
    CT_BUTTERFLY(aj5, aj13, zeta, t);
    CT_BUTTERFLY(aj6, aj14, zeta, t);
    CT_BUTTERFLY(aj7, aj15, zeta, t);
    // Second layer (len = 4), one root per half
    // First Half
    zeta = roots[k2];
    k2 += 1;
    CT_BUTTERFLY(aj0, aj4, zeta, t);
    CT_BUTTERFLY(aj1, aj5, zeta, t);
    CT_BUTTERFLY(aj2, aj6, zeta, t);
    CT_BUTTERFLY(aj3, aj7, zeta, t);
    // Second Half
    zeta = roots[k2];
    k2 += 1;
    CT_BUTTERFLY(aj8, aj12, zeta, t);
    CT_BUTTERFLY(aj9, aj13, zeta, t);
    CT_BUTTERFLY(aj10, aj14, zeta, t);
    CT_BUTTERFLY(aj11, aj15, zeta, t);
    // Third layer (len = 2), one root per quarter
    // First Quarter
    zeta = roots[k3];
    k3 += 1;
    CT_BUTTERFLY(aj0, aj2, zeta, t);
    CT_BUTTERFLY(aj1, aj3, zeta, t);
    // Second Quarter
    zeta = roots[k3];
    k3 += 1;
    CT_BUTTERFLY(aj4, aj6, zeta, t);
    CT_BUTTERFLY(aj5, aj7, zeta, t);
    // Third Quarter
    zeta = roots[k3];
    k3 += 1;
    CT_BUTTERFLY(aj8, aj10, zeta, t);
    CT_BUTTERFLY(aj9, aj11, zeta, t);
    // Forth Quarter
    zeta = roots[k3];
    k3 += 1;
    CT_BUTTERFLY(aj12, aj14, zeta, t);
    CT_BUTTERFLY(aj13, aj15, zeta, t);
    // Forth layer (len = 1), one root per eigth
    // First Eigth
    zeta = roots[k4];
    k4 += 1;
    CT_BUTTERFLY(aj0, aj1, zeta, t);
    // Second Eigth
    zeta = roots[k4];
    k4 += 1;
    CT_BUTTERFLY(aj2, aj3, zeta, t);
    // Third Eigth
    zeta = roots[k4];
    k4 += 1;
    CT_BUTTERFLY(aj4, aj5, zeta, t);
    // Forth Eigth
    zeta = roots[k4];
    k4 += 1;
    CT_BUTTERFLY(aj6, aj7, zeta, t);
    // Fifth Eigth
    zeta = roots[k4];
    k4 += 1;
    CT_BUTTERFLY(aj8, aj9, zeta, t);
    // Sixth Eigth
    zeta = roots[k4];
    k4 += 1;
    CT_BUTTERFLY(aj10, aj11, zeta, t);
    // Seventh Eigth
    zeta = roots[k4];
    k4 += 1;
    CT_BUTTERFLY(aj12, aj13, zeta, t);
    // Eigth Eigth
    zeta = roots[k4];
    k4 += 1;
    CT_BUTTERFLY(aj14, aj15, zeta, t);
    // Store
    a[j + 0] = aj0;
    a[j + 1] = aj1;
    a[j + 2] = aj2;
    a[j + 3] = aj3;
    a[j + 4] = aj4;
    a[j + 5] = aj5;
    a[j + 6] = aj6;
    a[j + 7] = aj7;
    a[j + 8] = aj8;
    a[j + 9] = aj9;
    a[j + 10] = aj10;
    a[j + 11] = aj11;
    a[j + 12] = aj12;
    a[j + 13] = aj13;
    a[j + 14] = aj14;
    a[j + 15] = aj15;
  }
}

#define GS_BUTTERFLY(aj, ak, zeta, t)                 \
  do {                                                \
    (t) = (aj) - (ak);                                \
    (aj) = mod_ct((aj) + (ak), (int64_t)q, qRq);      \
    (ak) = mod_mul_ct((zeta), (t), (int64_t)q, qinv); \
  } while (0)

void poly_intt(poly_t a, const poly_t roots, uint64_t q, int64_t qRq,
               int64_t qinv, uint32_t n) {
  int64_t t, zeta;
  // Merge first 4 layers using 16+ registers
  uint32_t k1 = n / 2;
  uint32_t k2 = n / 4;
  uint32_t k3 = n / 8;
  uint32_t k4 = n / 16;
  for (uint32_t j = 0; j < n; j += 16) {
    // Load
    int64_t aj0 = a[j + 0];
    int64_t aj1 = a[j + 1];
    int64_t aj2 = a[j + 2];
    int64_t aj3 = a[j + 3];
    int64_t aj4 = a[j + 4];
    int64_t aj5 = a[j + 5];
    int64_t aj6 = a[j + 6];
    int64_t aj7 = a[j + 7];
    int64_t aj8 = a[j + 8];
    int64_t aj9 = a[j + 9];
    int64_t aj10 = a[j + 10];
    int64_t aj11 = a[j + 11];
    int64_t aj12 = a[j + 12];
    int64_t aj13 = a[j + 13];
    int64_t aj14 = a[j + 14];
    int64_t aj15 = a[j + 15];
    // First layer (len = 1), one root per eigth
    // First Eigth
    zeta = roots[k1];
    k1 += 1;
    GS_BUTTERFLY(aj0, aj1, zeta, t);
    // Second Eigth
    zeta = roots[k1];
    k1 += 1;
    GS_BUTTERFLY(aj2, aj3, zeta, t);
    // Third Eigth
    zeta = roots[k1];
    k1 += 1;
    GS_BUTTERFLY(aj4, aj5, zeta, t);
    // Forth Eigth
    zeta = roots[k1];
    k1 += 1;
    GS_BUTTERFLY(aj6, aj7, zeta, t);
    // Fifth Eigth
    zeta = roots[k1];
    k1 += 1;
    GS_BUTTERFLY(aj8, aj9, zeta, t);
    // Sixth Eigth
    zeta = roots[k1];
    k1 += 1;
    GS_BUTTERFLY(aj10, aj11, zeta, t);
    // Seventh Eigth
    zeta = roots[k1];
    k1 += 1;
    GS_BUTTERFLY(aj12, aj13, zeta, t);
    // Eigth Eigth
    zeta = roots[k1];
    k1 += 1;
    GS_BUTTERFLY(aj14, aj15, zeta, t);
    // Second layer (len = 2), one root per quarter
    // First Quarter
    zeta = roots[k2];
    k2 += 1;
    GS_BUTTERFLY(aj0, aj2, zeta, t);
    GS_BUTTERFLY(aj1, aj3, zeta, t);
    // Second Quarter
    zeta = roots[k2];
    k2 += 1;
    GS_BUTTERFLY(aj4, aj6, zeta, t);
    GS_BUTTERFLY(aj5, aj7, zeta, t);
    // Third Quarter
    zeta = roots[k2];
    k2 += 1;
    GS_BUTTERFLY(aj8, aj10, zeta, t);
    GS_BUTTERFLY(aj9, aj11, zeta, t);
    // Forth Quarter
    zeta = roots[k2];
    k2 += 1;
    GS_BUTTERFLY(aj12, aj14, zeta, t);
    GS_BUTTERFLY(aj13, aj15, zeta, t);
    // Third layer (len = 4), one root per half
    // First Half
    zeta = roots[k3];
    k3 += 1;
    GS_BUTTERFLY(aj0, aj4, zeta, t);
    GS_BUTTERFLY(aj1, aj5, zeta, t);
    GS_BUTTERFLY(aj2, aj6, zeta, t);
    GS_BUTTERFLY(aj3, aj7, zeta, t);
    // Second Half
    zeta = roots[k3];
    k3 += 1;
    GS_BUTTERFLY(aj8, aj12, zeta, t);
    GS_BUTTERFLY(aj9, aj13, zeta, t);
    GS_BUTTERFLY(aj10, aj14, zeta, t);
    GS_BUTTERFLY(aj11, aj15, zeta, t);
    // Forth layer (len = 8), one root
    zeta = roots[k4];
    k4 += 1;
    GS_BUTTERFLY(aj0, aj8, zeta, t);
    GS_BUTTERFLY(aj1, aj9, zeta, t);
    GS_BUTTERFLY(aj2, aj10, zeta, t);
    GS_BUTTERFLY(aj3, aj11, zeta, t);
    GS_BUTTERFLY(aj4, aj12, zeta, t);
    GS_BUTTERFLY(aj5, aj13, zeta, t);
    GS_BUTTERFLY(aj6, aj14, zeta, t);
    GS_BUTTERFLY(aj7, aj15, zeta, t);
    // Store
    a[j + 0] = aj0;
    a[j + 1] = aj1;
    a[j + 2] = aj2;
    a[j + 3] = aj3;
    a[j + 4] = aj4;
    a[j + 5] = aj5;
    a[j + 6] = aj6;
    a[j + 7] = aj7;
    a[j + 8] = aj8;
    a[j + 9] = aj9;
    a[j + 10] = aj10;
    a[j + 11] = aj11;
    a[j + 12] = aj12;
    a[j + 13] = aj13;
    a[j + 14] = aj14;
    a[j + 15] = aj15;
  }

  uint32_t log2_n = floor_log2(n / 16);
  uint32_t n3_rem = log2_n % 3;
  uint32_t n3 = 1 << (log2_n + 4 - n3_rem);
  // Merge 3 layers using 16+ registers
  for (uint32_t len = 16; len < n3; len <<= 3) {
    k1 = n / (2 * len);
    k2 = k1 >> 1;
    k3 = k1 >> 2;
    for (uint32_t start = 0; start < n; start += 8 * len) {
      // Roots
      int64_t w10 = roots[k1];
      int64_t w11 = roots[k1 + 1];
      int64_t w12 = roots[k1 + 2];
      int64_t w13 = roots[k1 + 3];
      k1 += 4;
      int64_t w20 = roots[k2];
      int64_t w21 = roots[k2 + 1];
      k2 += 2;
      int64_t w3 = roots[k3];
      k3 += 1;
      for (uint32_t j = start; j < start + len; j++) {
        // Load
        int64_t aj0 = a[j + 0 * len];
        int64_t aj1 = a[j + 1 * len];
        int64_t aj2 = a[j + 2 * len];
        int64_t aj3 = a[j + 3 * len];
        int64_t aj4 = a[j + 4 * len];
        int64_t aj5 = a[j + 5 * len];
        int64_t aj6 = a[j + 6 * len];
        int64_t aj7 = a[j + 7 * len];
        // First layer (len), one root per quarter
        GS_BUTTERFLY(aj0, aj1, w10, t);
        GS_BUTTERFLY(aj2, aj3, w11, t);
        GS_BUTTERFLY(aj4, aj5, w12, t);
        GS_BUTTERFLY(aj6, aj7, w13, t);
        // Second layer (2*len), one root per half
        GS_BUTTERFLY(aj0, aj2, w20, t);
        GS_BUTTERFLY(aj1, aj3, w20, t);
        GS_BUTTERFLY(aj4, aj6, w21, t);
        GS_BUTTERFLY(aj5, aj7, w21, t);
        // Third layer (4*len), one root
        GS_BUTTERFLY(aj0, aj4, w3, t);
        GS_BUTTERFLY(aj1, aj5, w3, t);
        GS_BUTTERFLY(aj2, aj6, w3, t);
        GS_BUTTERFLY(aj3, aj7, w3, t);
        // Store
        a[j + 0 * len] = aj0;
        a[j + 1 * len] = aj1;
        a[j + 2 * len] = aj2;
        a[j + 3 * len] = aj3;
        a[j + 4 * len] = aj4;
        a[j + 5 * len] = aj5;
        a[j + 6 * len] = aj6;
        a[j + 7 * len] = aj7;
      }
    }
  }

  if (n3_rem == 1) {  // Just last layer
    // Roots
    int64_t w = roots[1];
    for (uint32_t j = 0; j < n / 2; j++) {
      // Load
      int64_t aj0 = a[j + 0 * (n / 2)];
      int64_t aj1 = a[j + 1 * (n / 2)];
      // First Layer
      GS_BUTTERFLY(aj0, aj1, w, t);
      // Store
      a[j + 0 * (n / 2)] = aj0;
      a[j + 1 * (n / 2)] = aj1;
    }
  } else if (n3_rem == 2) {  // Two last layers
    // Roots
    int64_t w10 = roots[2];
    int64_t w11 = roots[3];
    int64_t w2 = roots[1];
    for (uint32_t j = 0; j < n / 4; j++) {
      // Load
      int64_t aj0 = a[j + 0 * (n / 4)];
      int64_t aj1 = a[j + 1 * (n / 4)];
      int64_t aj2 = a[j + 2 * (n / 4)];
      int64_t aj3 = a[j + 3 * (n / 4)];
      // First layer (n/4)
      GS_BUTTERFLY(aj0, aj1, w10, t);
      GS_BUTTERFLY(aj2, aj3, w11, t);
      // Second layer (n/2)
      GS_BUTTERFLY(aj0, aj2, w2, t);
      GS_BUTTERFLY(aj1, aj3, w2, t);
      // Store
      a[j + 0 * (n / 4)] = aj0;
      a[j + 1 * (n / 4)] = aj1;
      a[j + 2 * (n / 4)] = aj2;
      a[j + 3 * (n / 4)] = aj3;
    }
  }

  // R^2/n multiplication performed in poly_add/poly_rsub
}

void poly_gen_roots(poly_t roots, uint32_t n, int64_t q, int64_t qR,
                    int64_t qinv) {
  // Find one primitive root (n > 2)
  uint32_t log2_n = floor_log2(n);
  int64_t root0 = find_smallest_2nth_root(n, q, qR, qinv);
  root0 = mod_mul(root0, qR, q);  // montgomery factor
  roots[0] = 1;
  int64_t etan = root0;
  roots[n / 2] = etan;
  for (uint32_t i = 2; i < n; i++) {
    etan = mod_mul_ct(etan, root0, q, qinv);
    roots[bit_reverse32(i, log2_n)] = etan;
  }
}

void poly_gen_inverse_roots(poly_t roots, uint32_t n, int64_t q, int64_t qinv) {
  // Find one primitive root (n > 2)
  uint32_t log2_n = floor_log2(n);
  int64_t etan = -1;             // x^n
  int64_t root0 = roots[n / 2];  // Previously poly_gen_roots output
  roots[0] = 1;
  etan = mod_mul(etan, root0, q);
  roots[n - 1] = etan;  // bit_reverse32(n - 1, log2_n)
  for (uint32_t i = 2; i < n; i++) {
    etan = mod_mul_ct(etan, root0, q, qinv);
    roots[bit_reverse32(n - i, log2_n)] = etan;
  }
}
