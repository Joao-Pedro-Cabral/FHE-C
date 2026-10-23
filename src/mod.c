
#include "mod.h"

#include <stdio.h>
#include <stdlib.h>

#include "utils.h"

int64_t mod_pow_ct(int64_t a, uint64_t b, int64_t q, int64_t qR, int64_t qinv) {
  int64_t c = 1, d = mod_mul(a, qR, q);
  while (b > 0) {
    if ((b & 0x1) == 1) c = mod_mul_ct(c, d, q, qinv);
    d = mod_mul_ct(d, d, q, qinv);
    b >>= 1;
  }
  return c;
}

uint64_t mod64_pow_ct(uint64_t a, uint64_t b) {
  uint64_t c = 1, d = a;
  while (b > 0) {
    if ((b & 0x1) == 1) c = c * d;
    d = d * d;
    b >>= 1;
  }
  return c;
}

static uint64_t mod_pow0(uint64_t a, uint64_t b, uint64_t q) {
  uint64_t c = 1, d = a;
  while (b > 0) {
    if ((b & 0x1) == 1) c = mod_mul0(c, d, q);
    d = mod_mul0(d, d, q);
    b >>= 1;
  }
  return c;
}

uint32_t mod_pow2(uint32_t a, uint32_t b, uint32_t q) {
  uint32_t c = 1, d = a;
  uint32_t q_1 = q - 1;
  while (b > 0) {
    if ((b & 0x1) == 1) c = (c * d) & q_1;
    d = (d * d) & q_1;
    b >>= 1;
  }
  return c;
}

static inline uint64_t gcd(uint64_t a, uint64_t b) {
  while (b != 0) {
    uint64_t r = a % b;
    a = b;
    b = r;
  }
  return a;
}

#define MILLER_N_TESTS 7
#define MILLER_BASIS_PRIMES 6

// Enough for p < 2^64
static const uint64_t miller_basis[MILLER_N_TESTS] = {
    2UL, 325UL, 9375UL, 28178UL, 450775UL, 9780504UL, 1795265022UL};

// We need to check these guys first to prevent these primes being considered as
// composite
static const uint64_t miller_base_primes[MILLER_BASIS_PRIMES] = {
    13UL, 19UL, 73UL, 193UL, 407521UL, 299210837UL};

// We suppose that the input is on the form d2^s+1
static bool miller_rabin_primality_test(uint64_t p) {
  uint32_t s = (uint32_t)__builtin_ctzll(p - 1);
  uint64_t d = (p - 1) >> s;
  // Check some base cases
  if ((p % 2 == 0) || (p % 3 == 0) || (p % 5 == 0)) return false;
  for (uint32_t i = 0; i < MILLER_BASIS_PRIMES; i++) {
    if (p == miller_base_primes[i]) return true;
  }
  // Test
  for (uint32_t i = 0; i < MILLER_N_TESTS; i++) {
    uint64_t a = miller_basis[i];
    uint64_t x = mod_pow0(a, d, p);
    bool composity;
    for (uint32_t j = 0; j < s; j++) {
      // p-1 and -p+1 are impossible because 4*p < 2^64
      composity = (x != 1 && x != p - 1);
      x = mod_mul0(x, x, p);
      if (x == 1 && composity) return false;
    }
    if (x != 1) return false;
  }
  return true;
}

static inline uint64_t pollard_func(uint64_t y, uint64_t a, uint64_t p) {
  // y = y^2 + a mod p
  y = mod_mul0(y, y, p);
  y += a;
  if (y > p) y -= p;
  return y;
}

static uint64_t brent_pollard_factor(uint64_t p) {
  const uint64_t m = 1000;
  uint64_t a, x, y, ys, r, q, g;
  a = 1;
  y = 1;
  r = 1;
  q = 1;

  do {
    x = y;
    for (uint64_t i = 0; i < r; i++) {
      y = pollard_func(y, a, p);
    }

    uint64_t k = 0;
    do {
      for (uint64_t i = 0; i < m && i < r - k; i++) {
        ys = y;
        y = pollard_func(y, a, p);
        // q = q |x-y| mod p
        q = mod_mul0(q, (x > y) ? x - y : y - x, p);
      }
      g = gcd(q, p);
      k += m;
    } while (k < r && g == 1);

    r <<= 1;
  } while (g == 1);

  if (g == p) {
    do {
      ys = pollard_func(ys, a, p);
      g = gcd((x > ys) ? x - ys : ys - x, p);
    } while (g == 1);
  }

  return g;
}

static uint32_t find_primes_factors_rec(uint64_t* s, uint32_t k, uint64_t m) {
  if (m == 1) {
    return k;
  }
  if (miller_rabin_primality_test(m)) {
    s[k] = m;
    k++;
    return k;
  }

  uint64_t factor = brent_pollard_factor(m);
  k = find_primes_factors_rec(s, k, factor);
  k = find_primes_factors_rec(s, k, m / factor);
  return k;
}

#define BASE_FACTORS_TESTS 7

static const uint32_t prime_factors[BASE_FACTORS_TESTS] = {3,  5,  7, 11,
                                                           13, 17, 19};

static uint32_t find_primes_factors(uint64_t* s, uint64_t m) {
  if (m == 0 || m == 1) return 0;
  // Remove 2 factors
  uint32_t k = 0;
  if ((m & 0x1) == 0) {
    uint32_t m2 = (uint32_t)__builtin_ctzll(m);
    m = m >> m2;
    s[0] = 2;
    k = 1;
  }
  // Remove some other base cases
  for (uint64_t i = 0; i < BASE_FACTORS_TESTS; i++) {
    uint32_t prime = (uint64_t)prime_factors[i];
    if (m % prime == 0) {
      s[k] = prime;
      k++;
      do {
        m /= prime;
      } while (m % prime == 0);
    }
  }
  // Find larger factors
  return find_primes_factors_rec(s, k, m);
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
      if (mod_pow_ct(r, primes[j], q, qR, qinv) == 1) {
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
      mod_pow_ct((int64_t)q_1th_root, ((uint64_t)q - 1) / (2 * n), q, qR, qinv);
  // To [0, 2q[
  if (smallest < 0) smallest += q;
  int64_t root2 = (int64_t)(((__int128_t)smallest * qR) % q);  // xR
  root2 = mod_mul_ct(root2, root2, q, qinv);                   // x^2R
  int64_t root = smallest;
  // Find smallest nth root (OpenFHE API)
  for (uint32_t i = 1; i < n; i++) {
    // mod_mul_ct is centered, but here I want [0, q[
    root = mod_mul_ct(root, root2, q, qinv);  // xy
    if (root < 0) root += q;
    if ((root < smallest) && (root != 1)) smallest = root;
  }
  // Centered
  if (smallest > q / 2) smallest -= q;
  return smallest;
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
