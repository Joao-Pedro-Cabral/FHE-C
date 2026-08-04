
#ifndef __UTILS_H__
#define __UTILS_H__

static inline int ceil_log2(int n) {
  if (n == 1) return 0;

  return 32 - __builtin_clz(n - 1);
}

static inline int floor_log2(int n) {
  return 31 - __builtin_clz(n);
}

static inline int bit_reverse32(int x) {
  return __builtin_bitreverse32(x);
}

#endif
