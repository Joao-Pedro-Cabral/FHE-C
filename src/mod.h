
#ifndef __MOD_H__
#define __MOD_H__

#include <stdbool.h>
#include <stdint.h>

// We do static inline here for remove function call overhead

// Modular Inverse
static inline int64_t mod_inverse(int64_t a, int64_t q) {
  if (a < 0) a += q;  // Function for [0, ..., q-1]
  int64_t t = 0, new_t = 1;
  int64_t r = q, new_r = a % q;

  while (new_r != 0) {
    int64_t quotient = r / new_r;

    int64_t tmp = new_t;
    new_t = t - quotient * new_t;
    t = tmp;

    tmp = new_r;
    new_r = r - quotient * new_r;
    r = tmp;
  }

  if (t > q / 2) t -= q;
  if (t < -q / 2) t += q;

  return t;
}

// Modular Reduction
static inline int64_t mod(int64_t a, int64_t q) {
  a = a % q;
  if (a > q / 2) a -= q;
  if (a < -q / 2) a += q;
  return a;
}

// Modular Multiplication
static inline int64_t mod_mul(int64_t a, int64_t b, int64_t q) {
  __int128_t prod = (__int128_t)a * b;
  a = (int64_t)(prod % q);
  if (a > q / 2) a -= q;
  if (a < -q / 2) a += q;
  return a;
}

#endif  // __MOD_H__
