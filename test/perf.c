
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "ckks.h"
#include "poly.h"
#include "random.h"
#include "rlwe.h"
#include "vector.h"
#include "xof.h"

#ifndef N2
#define N2 16
#endif

#ifndef L
#define L 4
#endif

#ifndef N
#define N 512
#endif

#ifndef N_TESTS
#define N_TESTS 1000
#endif

// "Don't care" for performance
static const uint64_t q0_bits = 60;
static const uint64_t scale_bits = 30;

static uint64_t time_ns(void) {
  struct timespec ts;
  clock_gettime(CLOCK_MONOTONIC, &ts);

  return (uint64_t)ts.tv_sec * 1000000000ULL + (uint64_t)ts.tv_nsec;
}

static void perf_elem_wise(vec_t* v, poly_t d, const poly_t a, const poly_t b,
                           const poly_t c, int64_t x, uint32_t h, uint64_t q,
                           int64_t R_q, int64_t qinvmR, uint64_t scale,
                           uint32_t n2, uint32_t n) {
  uint64_t start = time_ns();

#pragma omp parallel
  {
    for (uint32_t i = 0; i < N_TESTS; i++) {
      vector_to_poly(d, v, scale, n2, n);
      poly_to_vector(v, d, scale, n2, n);
      poly_pointwise_macc_exact(d, a, b, c, q, qinvmR, R_q, n);
      poly_permute_omp(d, a, h, n);
      poly_pointwise_mul_omp(d, a, b, q, qinvmR, n);
      poly_mulx_omp(d, x, q, qinvmR, n);
    }
  }

  uint64_t end = time_ns();
  double total = (double)(end - start) / 1000.0;

  printf("Element-Wise) Total: %.3f us, AVG: %.3f us\n", total,
         total / N_TESTS);
}

static void perf_ternary_sampler(poly_t s, xof_state_t* xof_state,
                                 uint64_t domain_separator, uint8_t* buf,
                                 uint32_t n) {
  uint64_t start = time_ns();

#pragma omp parallel
  {
    for (uint32_t i = 0; i < N_TESTS; i++) {
      poly_ternary_sampler(s, xof_state, domain_separator, buf, n);
    }
  }

  uint64_t end = time_ns();
  double total = (double)(end - start) / 1000.0;

  printf("Ternary Sampler) Total: %.3f us, AVG: %.3f us\n", total,
         total / N_TESTS);
}

static void perf_cbd_add(poly_t p, uint8_t* buf, const poly_t m,
                         xof_state_t* xof) {
  uint64_t start = time_ns();

#pragma omp parallel
  {
    for (uint32_t i = 0; i < N_TESTS; i++) {
      poly_cbd24_sampler_add(p, xof, 0, buf, m, N);
    }
  }

  uint64_t end = time_ns();

  double total = (double)(end - start) / 1000.0;

  printf("CBD) Total: %.3f us, AVG: %.3f us\n", total, total / N_TESTS);
}

static void perf_fft(vec_t* v, const vec_t* w) {
  uint64_t start = time_ns();

#pragma omp parallel
  {
    for (uint32_t i = 0; i < N_TESTS; i++) {
      fft(v, w, N2);
    }
  }

  uint64_t end = time_ns();

  double total = (double)(end - start) / 1000.0;

  printf("FFT) Total: %.3f us, AVG: %.3f us\n", total, total / N_TESTS);
}

static void perf_ifft(vec_t* v, const vec_t* w) {
  uint64_t start = time_ns();

#pragma omp parallel
  {
    for (uint32_t i = 0; i < N_TESTS; i++) {
      ifft(v, w, N2);
    }
  }

  uint64_t end = time_ns();

  double total = (double)(end - start) / 1000.0;

  printf("IFFT) Total: %.3f us, AVG: %.3f us\n", total, total / N_TESTS);
}

static void perf_ntt(poly_t p, const poly_t r, uint64_t q, int64_t R_q,
                     int64_t qinvmR) {
  uint64_t start = time_ns();

  for (uint32_t i = 0; i < N_TESTS; i++) {
    poly_ntt(p, r, q, R_q, qinvmR, N);
  }

  uint64_t end = time_ns();

  double total = (double)(end - start) / 1000.0;

  printf("NTT) Total: %.3f us, AVG: %.3f us\n", total, total / N_TESTS);
}

static void perf_intt(poly_t p, const poly_t r, uint64_t q, int64_t R_q,
                      int64_t qinvmR) {
  for (uint32_t i = 0; i < N; i++) {
    p[i] = i;
  }

  uint64_t start = time_ns();

#pragma omp parallel
  {
    for (uint32_t i = 0; i < N_TESTS; i++) {
      poly_intt(p, r, q, R_q, qinvmR, N);
    }
  }

  uint64_t end = time_ns();

  double total = (double)(end - start) / 1000.0;

  printf("INTT) Total: %.3f us, AVG: %.3f us\n", total, total / N_TESTS);
}

static void perf_context_generation(ckks_params_t* params) {
  uint64_t start = time_ns();

  for (uint32_t i = 0; i < N_TESTS; i++) {
    ckks_context_generation(params, q0_bits, scale_bits, N2, L, N, 1, NULL, 0);
  }

  uint64_t end = time_ns();

  double total = (double)(end - start) / 1000.0;

  printf("ContextGen) Total: %.3f us, AVG: %.3f us\n", total, total / N_TESTS);
}

static void perf_keygen(ckks_params_t* params) {
  uint64_t start = time_ns();

  for (uint32_t i = 0; i < N_TESTS; i++) {
    ckks_secret_keygen(params);
    ckks_public_keygen(params);
  }

  uint64_t end = time_ns();

  double total = (double)(end - start) / 1000.0;

  printf("KeyGen) Total: %.3f us, AVG: %.3f us\n", total, total / N_TESTS);
}

static void perf_relinearization_keygen(ckks_params_t* params) {
  uint64_t start = time_ns();

  for (uint32_t i = 0; i < N_TESTS; i++) {
    ckks_relinearization_keygen(params);
  }

  uint64_t end = time_ns();

  double total = (double)(end - start) / 1000.0;

  printf("Relinearization KeyGen) Total: %.3f us, AVG: %.3f us\n", total,
         total / N_TESTS);
}

static void perf_conjugation_keygen(ckks_params_t* params) {
  uint64_t start = time_ns();

  for (uint32_t i = 0; i < N_TESTS; i++) {
    ckks_conjugation_keygen(params);
  }

  uint64_t end = time_ns();

  double total = (double)(end - start) / 1000.0;

  printf("Conjugation KeyGen) Total: %.3f us, AVG: %.3f us\n", total,
         total / N_TESTS);
}

static void perf_inner_product_keygen(ckks_params_t* params) {
  uint64_t start = time_ns();

  for (uint32_t i = 0; i < N_TESTS; i++) {
    ckks_inner_product_keygen(params);
  }

  uint64_t end = time_ns();

  double total = (double)(end - start) / 1000.0;

  printf("Inner Product KeyGen) Total: %.3f us, AVG: %.3f us\n", total,
         total / N_TESTS);
}

static void perf_encode(poly_t m, const vec_t* v, const ckks_params_t* params) {
  uint64_t start = time_ns();

  for (uint32_t i = 0; i < N_TESTS; i++) {
    ckks_encode(m, v, params);
  }

  uint64_t end = time_ns();

  double total = (double)(end - start) / 1000.0;

  printf("Encode) Total: %.3f us, AVG: %.3f us\n", total, total / N_TESTS);
}

static void perf_decode(vec_t* v, const poly_t m, uint32_t level,
                        const ckks_params_t* params) {
  uint64_t start = time_ns();

  for (uint32_t i = 0; i < N_TESTS; i++) {
    ckks_decode(v, m, level, params);
  }

  uint64_t end = time_ns();

  double total = (double)(end - start) / 1000.0;

  printf("Decode) Total: %.3f us, AVG: %.3f us\n", total, total / N_TESTS);
}

static void perf_encrypt(ciphertext_t* ct, poly_t m, ckks_params_t* params) {
  uint64_t start = time_ns();

  for (uint32_t i = 0; i < N_TESTS; i++) {
    ckks_encrypt(ct, m, params);
  }

  uint64_t end = time_ns();

  double total = (double)(end - start) / 1000.0;

  printf("Encrypt) Total: %.3f us, AVG: %.3f us\n", total, total / N_TESTS);
}

static void perf_decrypt(poly_t m, const ciphertext_t* ct,
                         const ckks_params_t* params) {
  uint64_t start = time_ns();

  for (uint32_t i = 0; i < N_TESTS; i++) {
    ckks_decrypt(m, ct, params);
  }

  uint64_t end = time_ns();

  double total = (double)(end - start) / 1000.0;

  printf("Decrypt) Total: %.3f us, AVG: %.3f us\n", total, total / N_TESTS);
}

int main() {
  vec_t v, w;
  poly_t p, a, b, c, r, r_1, m;
  xof_state_t xof_state;
  uint8_t seed[SEED_BYTES];
  rns_params_t rns_params;
  ckks_params_t ckks_params;
  ciphertext_t ct;

  // Init
  randombytes(seed, SEED_BYTES);
  xof_absorb(&xof_state, seed);

  vec_new(&v, N2);
  vec_new(&w, N2);
  poly_new(&p, N);
  poly_new(&r, N);
  poly_new(&r_1, N);
  poly_new(&m, N);
  poly_new(&a, N);
  poly_new(&b, N);
  poly_new(&c, N);

  uint32_t buflen = (((4 * (N / 4) + 2) / 3) + XOF_RATE - 1) / XOF_RATE;
  buflen = (buflen > 6 * N) ? buflen : 6 * N;
  uint8_t* buf = malloc(buflen);

  rns_params_new(&rns_params, q0_bits, scale_bits, L, N, 1, NULL, 0);

  vec_gen_roots(&w, N2);
  poly_gen_roots(r, N, (int64_t)rns_params.q[0], rns_params.Rmq[0],
                 rns_params.qinvmR[0]);
  poly_gen_inverse_roots(r_1, r[N / 2], N, (int64_t)rns_params.q[0],
                         rns_params.qinvmR[0]);

  // Perf Main Operations
  perf_cbd_add(p, buf, m, &xof_state);
  perf_ternary_sampler(p, &xof_state, 0, buf, N);
  perf_elem_wise(&v, p, a, b, c, 10, 5, rns_params.q[0], rns_params.R_q[0],
                 rns_params.qinvmR[0], 1 << scale_bits, N2, N);
  perf_fft(&v, &w);
  perf_ifft(&v, &w);
  perf_ntt(p, r, rns_params.q[0], rns_params.R_q[0], rns_params.qinvmR[0]);
  perf_intt(p, r_1, rns_params.q[0], rns_params.R_q[0], rns_params.qinvmR[0]);

  // Perf CKKS Operations
  // Context generation
  perf_context_generation(&ckks_params);

  // Key generation operations
  perf_keygen(&ckks_params);
  perf_relinearization_keygen(&ckks_params);
  perf_conjugation_keygen(&ckks_params);
  perf_inner_product_keygen(&ckks_params);

  // Encoding and decoding
  perf_encode(p, &v, &ckks_params);
  perf_decode(&v, p, 1, &ckks_params);

  // Encryption and decryption
  ciphertext_new(&ct, &ckks_params.rns);
  perf_encrypt(&ct, p, &ckks_params);
  perf_decrypt(p, &ct, &ckks_params);

  // Free
  free(buf);
  vec_free(&v);
  vec_free(&w);
  poly_free(&p);
  poly_free(&r);
  poly_free(&r_1);
  poly_free(&a);
  poly_free(&b);
  poly_free(&c);
  poly_free(&m);
  rns_params_free(&rns_params);
  ckks_free(&ckks_params);
  ciphertext_free(&ct);
}
