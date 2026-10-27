
#ifndef __MOD_H__
#define __MOD_H__

#include <stdbool.h>
#include <stdint.h>

// We do static inline here for remove function call overhead

int64_t mod_pow_ct(int64_t a, uint64_t b, int64_t q, int64_t Rmq, int64_t qinvmR);

uint64_t mod64_pow_ct(uint64_t a, uint64_t b);

// Modular Multiplication (ct: constant-time)
static inline int64_t mod_mul_ct(int64_t a, int64_t b, int64_t q,
                                 int64_t qinvmR) {
  int64_t c0 = a * b;
  int64_t c1 = (int64_t)((__int128_t)a * b >> 64);
  int64_t t = c0 * qinvmR;
  t = (int64_t)((__int128_t)t * q >> 64);
  t = c1 - t;
  return t;
}

// Modular Reduction
static inline int64_t mod_ct(int64_t a, int64_t q, int64_t R_q) {
  int64_t c =
      (int64_t)(((__int128_t)a * R_q + (1UL << 63)) >> 64);  // Round even
  return a - c * q;
}

static inline int64_t q_to_q2(int64_t a, int64_t q) {
  int64_t diff1 = (q / 2) - a;
  int64_t diff2 = a - (-q / 2);

  int64_t mask1 = diff1 >> 63;
  int64_t mask2 = diff2 >> 63;

  // Apply modifications using bitwise AND
  a -= (mask1 & q);
  a += (mask2 & q);

  return a;
}

static inline int64_t mod_ct_exact(int64_t a, int64_t q, int64_t R_q) {
  int64_t c = mod_ct(a, q, R_q);

  return q_to_q2(c, q);
}

static inline int64_t mod_mul_ct_exact(int64_t a, int64_t b, int64_t q,
                                       int64_t qinvmR) {
  int64_t c = mod_mul_ct(a, b, q, qinvmR);

  return q_to_q2(c, q);
}

static inline uint64_t mod_mul0(uint64_t a, uint64_t b, uint64_t q) {
  __uint128_t prod = (__uint128_t)a * b;
  a = (uint64_t)(prod % q);
  return a;
}

static inline int64_t mod_mul(int64_t a, int64_t b, int64_t q) {
  __int128_t prod = (__int128_t)a * b;
  a = (int64_t)(prod % q);
  if (a > q / 2) a -= q;
  if (a < -q / 2) a += q;
  return a;
}

int64_t find_smallest_2nth_root(uint64_t n, int64_t q, int64_t Rmq,
                                int64_t qinvmR);

void find_primes(uint64_t* q, uint32_t q0_bits, uint32_t scale_bits, uint32_t n,
                 uint32_t l);

uint32_t mod_pow2(uint32_t a, uint32_t b, uint32_t q);

#endif  // __MOD_H__
