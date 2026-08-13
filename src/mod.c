
#include "mod.h"

#include <stdlib.h>

#include "utils.h"

int64_t mod_pow(int64_t a, uint64_t b, int64_t q, int64_t qR, int64_t qinv) {
  int64_t c = 1, d = mod_mul(a, qR, q);
  while (b > 0) {
    if ((b & 0x1) == 1) c = mod_mul_ct(c, d, q, qinv);
    d = mod_mul_ct(d, d, q, qinv);
    b >>= 1;
  }
  return c;
}

uint64_t mod64_pow(uint64_t a, uint64_t b) {
  uint64_t c = 1, d = a;
  while (b > 0) {
    if ((b & 0x1) == 1) c = c*d;
    d = d*d;
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

static uint64_t find_primitive_root(int64_t q, int64_t qR, int64_t qinv) {
  uint64_t log2_q =
      ceil64_log2((uint64_t)q);  // at maximum log2_q prime factors
  uint64_t* primes = malloc(log2_q * sizeof(uint64_t));
  int64_t phi = q - 1;
  int64_t r;
  uint32_t k = find_primes_factors(primes, (uint64_t)phi);
  for (uint32_t j = 0; j < k; j++) {
    primes[j] = (uint64_t)phi / primes[j];
  }
  for (r = 2; r <= phi; r++) {
    bool find = true;
    for (uint32_t j = 0; j < k; j++) {
      if (mod_pow(r, primes[j], q, qR, qinv) == 1) {
        find = false;
        break;
      }
    }
    if (find) break;
  }
  free(primes);
  return (uint64_t)r;
}

int64_t find_smallest_2nth_root(uint64_t n, int64_t q, int64_t qR, int64_t qinv) {
  uint64_t q_1th_root = find_primitive_root(q, qR, qinv);
  int64_t smallest =
      mod_pow((int64_t)q_1th_root, ((uint64_t)q - 1) / (2 * n), q, qR, qinv);
  // To [0, 2q[
  if (smallest < 0) smallest += q;
  int64_t root2 = (int64_t)(((__int128_t)smallest * qR) % q); // xR
  root2 = mod_mul_ct(root2, root2, q, qinv); // x^2R
  int64_t root = smallest;
  // Find smallest nth root (OpenFHE API)
  for (uint32_t i = 1; i < n; i++) {
    // mod_mul is centered, but here I want [0, q[
    root = mod_mul_ct(root, root2, q, qinv); // xy
    if (root < 0) root += q;
    if ((root < smallest) && (root != 1)) smallest = root;
  }
  // Centered
  if (smallest > q / 2) smallest -= q;
  return smallest;
}
