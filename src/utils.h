
#ifndef __UTILS_H__
#define __UTILS_H__

// We do static inline here for remove function call overhead

static inline uint64_t ceil64_log2(uint64_t n) {
  if (n == 1) return 0;

  return (uint64_t)(64 - __builtin_clzll(n - 1));
}

static inline uint32_t ceil_log2(uint32_t n) {
  if (n == 1) return 0;

  return (uint32_t)(32 - __builtin_clz(n - 1));
}

static inline uint32_t floor_log2(uint32_t n) {
  return (uint32_t)(31 - __builtin_clz(n));
}

static inline uint32_t bit_reverse32(uint32_t x, uint32_t k) {
  x = ((x >> 1) & 0x55555555u) | ((x & 0x55555555u) << 1);
  x = ((x >> 2) & 0x33333333u) | ((x & 0x33333333u) << 2);
  x = ((x >> 4) & 0x0F0F0F0Fu) | ((x & 0x0F0F0F0Fu) << 4);
  x = ((x >> 8) & 0x00FF00FFu) | ((x & 0x00FF00FFu) << 8);
  x = (x >> 16) | (x << 16);
  return x >> (32 - k);
}

#endif
