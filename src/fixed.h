
#ifndef __FIXED_H__
#define __FIXED_H__

#include <stdint.h>

// Fixed point mul
static inline int64_t fixed_mul(int64_t a, int64_t b) {
  __int128 p = (__int128)a * b;
  p += (INT64_C(1) << 62);
  return (int64_t)(p >> 63);  // exclude extra signal bit
}

// Fixed point add
static inline int64_t fixed_add(int64_t a, int64_t b) {
  __int128 p = (__int128)a + b;
  p += 1;
  return (int64_t)(p >> 1);
}

// Fixed point sub
static inline int64_t fixed_sub(int64_t a, int64_t b) {
  __int128 p = (__int128)a - b;
  p += 1;
  return (int64_t)(p >> 1);
}

// Fixed point shift right
static inline int64_t fixed_shift(int64_t a, int scale) {
  a += (INT64_C(1) << scale);
  return a;
}

#endif  // __FIXED_H__
