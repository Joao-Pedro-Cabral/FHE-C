
#include "mod.h"

#include <stdio.h>
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
    if ((b & 0x1) == 1) c = c * d;
    d = d * d;
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

int64_t find_smallest_2nth_root(uint64_t n, int64_t q, int64_t qR,
                                int64_t qinv) {
  uint64_t q_1th_root = find_primitive_root(q, qR, qinv);
  int64_t smallest =
      mod_pow((int64_t)q_1th_root, ((uint64_t)q - 1) / (2 * n), q, qR, qinv);
  // To [0, 2q[
  if (smallest < 0) smallest += q;
  int64_t root2 = (int64_t)(((__int128_t)smallest * qR) % q);  // xR
  root2 = mod_mul_ct(root2, root2, q, qinv);                   // x^2R
  int64_t root = smallest;
  // Find smallest nth root (OpenFHE API)
  for (uint32_t i = 1; i < n; i++) {
    // mod_mul is centered, but here I want [0, q[
    root = mod_mul_ct(root, root2, q, qinv);  // xy
    if (root < 0) root += q;
    if ((root < smallest) && (root != 1)) smallest = root;
  }
  // Centered
  if (smallest > q / 2) smallest -= q;
  return smallest;
}

#define MILLER_N_TESTS 7

// Enough for p < 2^64
uint64_t miller_primes[MILLER_N_TESTS] = {
    2UL, 325UL, 9375UL, 28178UL, 450775UL, 9780504UL, 1795265022UL};

// We suppose that the input is on the form d2^s+1
static bool miller_rabin_primality_test(uint64_t p) {
  uint32_t s = (uint32_t)__builtin_ctzll(p - 1);
  uint64_t d = (p - 1) >> s;
  // Check some base cases
  if ((p % 2 == 0) || (p % 3 == 0) || (p % 5 == 0)) return false;
  // Auxiliar values for fast modular operations
  int64_t pinv = (int64_t)mod64_pow(p, (1UL << 63) - 1);
  uint64_t R = (uint64_t)(((__uint128_t)1 << 64) - p);
  int64_t pR = (int64_t)(R % p);
  if (pR > (int64_t)(p / 2)) pR -= (int64_t)p;
  // Test
  for (uint32_t i = 0; i < MILLER_N_TESTS; i++) {
    uint64_t a = miller_primes[i];
    int64_t x = mod_pow((int64_t)a, d, (int64_t)p, pR, pinv);
    bool composity;
    for (uint32_t j = 0; j < s; j++) {
      // p-1 and -p+1 are impossible because 4*p < 2^64
      composity = (x != 1 && x != -1);
      x = mod_mul(x, x, (int64_t)p);
      if (x == 1 && composity) return false;
    }
  }
  return true;
}

// p_bits: number of bits of p
static uint64_t first_prime(uint32_t p_bits, uint32_t n) {
  // Since n is a power of two smaller than p_bits
  // Our initial guess is p = 2^p_bits + 2^(p_bits-1) + 1
  // We do like this to reduce the probability of failure to find the others
  // primes with the same bit-width
  uint64_t p = (1UL << (p_bits - 1)) + (1UL << (p_bits - 2)) + 1;
  while (!miller_rabin_primality_test(p)) {
    p += n;
  }
  if (ceil64_log2(p) != p_bits) {
    printf("First p too large: p = %ld > %ld\n", p,
           ((1UL << (p_bits - 1)) - 1));
  }
  return p;
}

// p_bits: number of bits of p
static uint64_t next_prime(uint64_t p, uint32_t n) {
  p = p + n;  // next number 1 mod n
  while (!miller_rabin_primality_test(p)) {
    p += n;
  }
  return p;
}

// p_bits: number of bits of p
static uint64_t previous_prime(uint64_t p, uint32_t n) {
  p = p - n;  // previous number 1 mod n
  while (!miller_rabin_primality_test(p)) {
    p -= n;
  }
  return p;
}

void find_primes(uint64_t* q, uint32_t q0_bits, uint32_t scale_bits, uint32_t n,
                 uint32_t l) {
  // First find the l scale_bits primes
  uint64_t q1 = first_prime(scale_bits, n);
  q[l - 1] = q1;
  double q_approx = (double)q1;
  uint64_t q_max = q1;
  uint32_t cnt = 0;
  for (uint32_t i = l - 2; i > 0; i--, cnt++) {
    q_approx = (double)(q_approx * q_approx) / (double)q[i + 1];
    uint64_t q_int = (uint64_t)(q_approx + 0.5);
    q_int = q_int - (q_int & (n - 1)) + 1;
    if (cnt % 2 == 0) {
      bool repeated = true;
      while (repeated) {
        repeated = false;
        q_int = previous_prime(q_int, n);
        for (size_t j = i + 1; j < l; j++) {
          if (q_int == q[j]) {
            repeated = true;
            break;
          }
        }
      }
    } else {
      bool repeated = true;
      while (repeated) {
        repeated = false;
        q_int = next_prime(q_int, n);
        for (size_t j = i + 1; j < l; j++) {
          if (q_int == q[j]) {
            repeated = true;
            break;
          }
        }
      }
    }
    q[i] = q_int;
    if (ceil64_log2(q_int) != scale_bits) {
      printf("q[%d] with incorrect bit-width: %ld instead of %d\n", i,
             ceil64_log2(q_int), scale_bits);
    }
    if (q_int > q_max) q_max = q_int;
  }
  // Find q0
  if (q0_bits == scale_bits) {
    q[0] = next_prime(q_max, n);
  } else {
    q[0] = first_prime(q0_bits, n);
  }
}
