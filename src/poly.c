
#include "poly.h"

#include <stdlib.h>

#include "fixed.h"
#include "mod.h"
#include "random.h"
#include "shake.h"
#include "utils.h"

void poly_new(poly* m, uint32_t n) { *m = malloc(n * sizeof(int64_t)); }

void poly_free(poly* m) { free(*m); }

void vector_to_poly(poly m, const vec* v, int64_t scale, uint32_t n) {
  for (uint32_t i = 0; i < n; i++) {
    m[i] = (int64_t)(v->rel[i] * (double)scale);
  }
}

void poly_to_vector(vec* v, const poly m, int64_t scale, uint32_t n) {
  // To increase the precision, we first compute (1 << int_part)/scale
  // Instead of /scale and *(1 << int_part)
  double scale2 = 1.0 / (double)scale;
  for (uint32_t i = 0; i < n; i++) {
    v->rel[i] = (double)m[i] * scale2;
    v->img[i] = 0;
  }
}

void poly_ternary_sampler(poly s, uint32_t n) {
  // Initialize PRNG
  uint8_t seed[SEED_BYTES];
  randombytes(seed, SEED_BYTES);
  shake128ctx state;
  shake128_absorb(&state, seed, SEED_BYTES);
  // Sample
  size_t blocks = (((4 * (n / 4) + 2) / 3) + SHAKE128_RATE - 1) / SHAKE128_RATE;
  size_t buflen = blocks * SHAKE128_RATE;
  uint8_t* buf = malloc(buflen);  // Initial guess
  shake128_squeezeblocks(buf, blocks, &state);
  uint32_t i = 0, j = 0;
  while (i < n) {
    // We're probably almost at the end so just fetch one more block
    if (j == buflen) {
      blocks = 1;
      buflen = SHAKE128_RATE;
      shake128_squeezeblocks(buf, blocks, &state);
      j = 0;
    }
    uint8_t x = buf[j];
    for (int k = 0; k < 8 / 2 && i < n; k++) {
      uint8_t y = (x >> (2 * k)) & 0x3;
      if (y < 3) {
        int8_t z = (int8_t)(y & 0x1) - (int8_t)(y >> 1);
        s[i] = z;
        i++;
      }
    }
    j++;
  }
  free(buf);
}

void poly_uniform_sampler(poly a, const uint8_t seed[SEED_BYTES], uint64_t q,
                          uint32_t n) {
  // Initialize PRNG
  shake128ctx state;
  shake128_absorb(&state, seed, SEED_BYTES);
  // Sample
  uint64_t bound = UINT64_C(-1) / q;
  bound *= q;
  // Normally q < 2^50, then the rejection probability (2^63-1-kq)/(2^63-1) is
  // extremely small, so I just fetch n int64_t plus the SHAKE rate rounding
  size_t blocks = (n * sizeof(int64_t));
  blocks = (blocks + SHAKE128_RATE - 1) / SHAKE128_RATE;
  size_t buflen = blocks * SHAKE128_RATE;
  // The rate is multiple of 8 so it's fine
  uint64_t* buf = (uint64_t*)malloc(buflen);  // Initial guess
  shake128_squeezeblocks((uint8_t*)buf, blocks, &state);
  uint32_t i = 0, j = 0;
  while (i < n) {
    // We're probably almost at the end so just fetch one more block
    if (j == buflen) {
      blocks = 1;
      buflen = SHAKE128_RATE;
      shake128_squeezeblocks((uint8_t*)buf, blocks, &state);
      j = 0;
    }
    uint64_t x = buf[j];
    if (x < bound) {
      int64_t y = mod((int64_t)x, (int64_t)q);
      a[i] = y;
      i++;
    }
    j++;
  }
  free(buf);
}

void poly_cbd24_sampler(poly e, uint32_t n) {
  uint8_t* buf = malloc(6 * n * sizeof(uint8_t));
  uint8_t seed[SEED_BYTES];
  randombytes(seed, SEED_BYTES);
  shake128(buf, 6 * n, seed, SEED_BYTES);
  for (uint32_t i = 0; i < n; i++) {
    int b1 = __builtin_popcount(buf[6 * i]);
    b1 += __builtin_popcount(buf[6 * i + 1]);
    b1 += __builtin_popcount(buf[6 * i + 2]);
    int b2 = __builtin_popcount(buf[6 * i + 3]);
    b2 += __builtin_popcount(buf[6 * i + 4]);
    b2 += __builtin_popcount(buf[6 * i + 5]);
    e[i] = b1 - b2;
  }
  free(buf);
}

void poly_add_no_mod(poly a, const poly b, uint32_t n) {
  for (uint32_t i = 0; i < n; i++) {
    a[i] = a[i] + b[i];
  }
}

void poly_add(poly a, const poly b, uint64_t q, uint32_t n) {
  for (uint32_t i = 0; i < n; i++) {
    a[i] = mod(a[i] + b[i], (int64_t)q);
  }
}

void poly_rsub(poly a, const poly b, uint64_t q, uint32_t n) {
  for (uint32_t i = 0; i < n; i++) {
    a[i] = mod(b[i] - a[i], (int64_t)q);
  }
}

void poly_pointwise_mul(poly a, const poly b, uint64_t q, uint32_t n) {
  for (uint32_t i = 0; i < n; i++) {
    a[i] = mod_mul(a[i], b[i], (int64_t)q);
  }
}

void poly_ntt(poly a, const roots eta, uint64_t q, uint32_t n) {
  int64_t zeta, t;

  for (uint32_t len = n / 2; len > 0; len >>= 1) {
    uint32_t k = n / (2 * len);
    for (uint32_t start = 0; start < n; start += 2 * len) {
      zeta = eta[k];
      k++;
      for (uint32_t j = start; j < start + len; ++j) {
        t = mod_mul(zeta, a[j + len], (int64_t)q);
        a[j + len] = mod(a[j] - t, (int64_t)q);
        a[j] = mod(a[j] + t, (int64_t)q);
      }
    }
  }
}

void poly_intt(poly a, const roots eta, uint64_t q, uint32_t n) {
  int64_t t, zeta;
  int64_t f = mod_inverse(n, (int64_t)q);  // mont^2/256

  for (uint32_t len = 1; len < n; len <<= 1) {
    uint32_t k = n / (2 * len);
    for (uint32_t start = 0; start < n; start += 2 * len) {
      zeta = eta[k];
      k++;
      for (uint32_t j = start; j < start + len; ++j) {
        t = a[j];
        a[j] = mod(t + a[j + len], (int64_t)q);
        a[j + len] = mod(t - a[j + len], (int64_t)q);
        a[j + len] = mod_mul(zeta, a[j + len], (int64_t)q);
      }
    }
  }

  for (uint32_t j = 0; j < n; ++j) {
    a[j] = mod_mul(a[j], f, (int64_t)q);
  }
}
