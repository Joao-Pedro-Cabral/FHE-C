
#include "mod.h"

#include <stdlib.h>

#include "utils.h"

int64_t mod_inverse(uint64_t a, uint64_t q) {
  int64_t t = 0, new_t = 1;
  uint64_t r = q, new_r = a % q;

  while (new_r != 0) {
    uint64_t quotient = r / new_r;

    int64_t tmp = new_t;
    new_t = t - (int64_t)quotient * new_t;
    t = tmp;

    uint64_t tmp2 = new_r;
    new_r = r - quotient * new_r;
    r = tmp2;
  }

  return t;
}

int64_t mod_inverse_big(__uint128_t a, __uint128_t q) {
  __int128_t t = 0, new_t = 1;
  __uint128_t r = q, new_r = a % q;

  while (new_r != 0) {
    __uint128_t quotient = r / new_r;

    __int128_t tmp = new_t;
    new_t = t - (__int128_t)quotient * new_t;
    t = tmp;

    __uint128_t tmp2 = new_r;
    new_r = r - quotient * new_r;
    r = tmp2;
  }

  return (int64_t)t;
}

static uint64_t mod_pow(uint64_t a, uint64_t b, uint64_t q) {
  uint64_t c = 1, d = a;
  while (b > 0) {
    if ((b & 0x1) == 1) c = mod_mulu(c, d, q);
    d = mod_mulu(d, d, q);
    b >>= 1;
  }
  return c;
}

static uint32_t find_primes_factors(uint64_t* s, uint64_t m) {
  uint32_t k = 0;
  for (uint64_t i = 2; i * i < m; i++) {
    if (m % i == 0) {
      s[k] = i;
      k++;
      do {
        m /= i;
      } while (m % i == 0);
    }
  }
  if (m > 1) {
    s[k] = m;
    k++;
  }
  return k;
}

static uint64_t find_primitive_root(int64_t q) {
  uint64_t log2_q =
      ceil64_log2((uint64_t)q);  // at maximum log2_q prime factors
  uint64_t* primes = malloc(log2_q * sizeof(uint64_t));
  uint64_t phi = (uint64_t)q - 1;
  uint64_t r;
  uint32_t k = find_primes_factors(primes, phi);
  for (r = 2; r <= phi; r++) {
    bool find = true;
    for (uint32_t j = 0; j < k; j++) {
      if (mod_pow(r, phi / primes[j], (uint64_t)q) == 1) {
        find = false;
        break;
      }
    }
    if (find) break;
  }
  free(primes);
  return r;
}

int64_t find_smallest_nth_root(uint64_t n, int64_t q) {
  uint64_t q_1th_root = find_primitive_root(q);
  uint64_t nth_root =
      mod_pow(q_1th_root, ((uint64_t)q - 1) / (2 * n), (uint64_t)q);
  int64_t smallest = (int64_t)nth_root;
  uint64_t root_pow2 = (uint64_t)(((__int128_t)nth_root * nth_root) % q);
  uint64_t root = nth_root;
  // Find smallest nth root (OpenFHE API)
  for (uint32_t i = 1; i < 2 * n; i++) {
    // mod_mul is centered, but here I want [0, q[
    root = (uint64_t)(((__int128_t)root * root_pow2) % q);
    if ((root < (uint64_t)smallest) && (root != 1)) smallest = (int64_t)root;
  }
  // Centered
  if (smallest > q / 2) smallest -= q;
  return smallest;
}
