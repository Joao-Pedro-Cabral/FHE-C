
#ifndef __MOD_H__
#define __MOD_H__

#include <stdbool.h>
#include <stdint.h>

// We do static inline here for remove function call overhead

// Modular Inverse
int64_t mod_inverse(uint64_t a, uint64_t q);

// Modular Multiplication
static inline int64_t mod_mul(int64_t a, int64_t b, int64_t q, int64_t qinv) {
  int64_t c0 = a * b;
  int64_t c1 = (int64_t)((__int128_t)a * b >> 64);
  int64_t t = c1 * qinv;
  t = (int64_t)((__int128_t)t * q >> 64);
  t = c0 - t;
  return t;
}

// Modular Reduction
static inline int64_t mod(int64_t a, int64_t q, int64_t qR, int64_t qinv) {
  return mod_mul(a, qR, q, qinv);
}

static inline int64_t mod_mul_exact(int64_t a, int64_t b, int64_t q,
                                    int64_t qinv) {
  int64_t c = mod_mul(a, b, q, qinv);

  int64_t diff1 = (q / 2) - c;
  int64_t diff2 = c - (-q / 2);

  int64_t mask1 = diff1 >> 63;
  int64_t mask2 = diff2 >> 63;

  // Apply modifications using bitwise AND
  c -= (mask1 & q);
  c += (mask2 & q);

  return c;
}

static inline int64_t mod_exact(int64_t a, int64_t q, int64_t qR,
                                int64_t qinv) {
  return mod_mul_exact(a, qR, q, qinv);
}

#endif  // __MOD_H__
