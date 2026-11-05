
#include <stdio.h>
#include <time.h>

#include "ckks.h"
#include "poly.h"
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

static void perf_cbd(poly_t p, xof_state_t* xof) {
  uint64_t start = time_ns();

  for (uint32_t i = 0; i < N_TESTS; i++) {
    poly_cbd24_sampler(p, xof, 0, N);
  }

  uint64_t end = time_ns();

  double total = (double)(end - start) / 1000.0;

  printf("CBD) Total: %.3f ms, AVG: %.3f ms\n", total, total / N_TESTS);
}

static void perf_fft(vec_t* v, const vec_t* w) {
  uint64_t start = time_ns();

  for (uint32_t i = 0; i < N_TESTS; i++) {
    fft(v, w, N2);
  }

  uint64_t end = time_ns();

  double total = (double)(end - start) / 1000.0;

  printf("FFT) Total: %.3f ms, AVG: %.3f ms\n", total, total / N_TESTS);
}

static void perf_ifft(vec_t* v, const vec_t* w) {
  uint64_t start = time_ns();

  for (uint32_t i = 0; i < N_TESTS; i++) {
    ifft(v, w, N2);
  }

  uint64_t end = time_ns();

  double total = (double)(end - start) / 1000.0;

  printf("IFFT) Total: %.3f ms, AVG: %.3f ms\n", total, total / N_TESTS);
}

static void perf_ntt(poly_t p, const poly_t r, uint64_t q, int64_t R_q,
                     int64_t qinvmR) {
  uint64_t start = time_ns();

  for (uint32_t i = 0; i < N_TESTS; i++) {
    poly_ntt(p, r, q, R_q, qinvmR, N);
  }

  uint64_t end = time_ns();

  double total = (double)(end - start) / 1000.0;

  printf("NTT) Total: %.3f ms, AVG: %.3f ms\n", total, total / N_TESTS);
}

static void perf_intt(poly_t p, const poly_t r, uint64_t q, int64_t R_q,
                      int64_t qinvmR) {
  uint64_t start = time_ns();

  for (uint32_t i = 0; i < N_TESTS; i++) {
    poly_intt(p, r, q, R_q, qinvmR, N);
  }

  uint64_t end = time_ns();

  double total = (double)(end - start) / 1000.0;

  printf("INTT) Total: %.3f ms, AVG: %.3f ms\n", total, total / N_TESTS);
}

static void perf_keygen(ckks_params_t* params) {
  uint64_t start = time_ns();

  for (uint32_t i = 0; i < (N_TESTS + 9) / 10; i++) {
    ckks_keygen(params, q0_bits, scale_bits, N2, L, N, 1, NULL, 0);
  }

  uint64_t end = time_ns();

  double total = (double)(end - start) / 1000.0;

  printf("KeyGen) Total: %.3f ms, AVG: %.3f ms\n", total, total / N_TESTS);
}

static void perf_encrypt(ciphertext_t* ct, poly_t m, const vec_t* v,
                         ckks_params_t* params) {
  uint64_t start = time_ns();

  for (uint32_t i = 0; i < N_TESTS; i++) {
    ckks_encode(m, v, params);
    ckks_encrypt(ct, m, params);
  }

  uint64_t end = time_ns();

  double total = (double)(end - start) / 1000.0;

  printf("Encrypt) Total: %.3f ms, AVG: %.3f ms\n", total, total / N_TESTS);
}

static void perf_decrypt(vec_t* v, poly_t m, const ciphertext_t* ct,
                         const ckks_params_t* params) {
  uint64_t start = time_ns();

  for (uint32_t i = 0; i < N_TESTS; i++) {
    ckks_decrypt(m, ct, params);
    ckks_decode(v, m, ct->level, params);
  }

  uint64_t end = time_ns();

  double total = (double)(end - start) / 1000.0;

  printf("Decrypt) Total: %.3f ms, AVG: %.3f ms\n", total, total / N_TESTS);
}

int main() {
  vec_t v, w;
  poly_t p, r, r_1;
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

  rns_params_new(&rns_params, q0_bits, scale_bits, L, N, 1, NULL, 0);

  vec_gen_roots(&w, N2);
  poly_gen_roots(r, N, (int64_t)rns_params.q[0], rns_params.Rmq[0],
                 rns_params.qinvmR[0]);
  poly_gen_inverse_roots(r_1, r[N / 2], N, (int64_t)rns_params.q[0],
                         rns_params.qinvmR[0]);

  // Perf Main Operations
  perf_cbd(p, &xof_state);
  perf_fft(&v, &w);
  perf_ifft(&v, &w);
  perf_ntt(p, r, rns_params.q[0], rns_params.R_q[0], rns_params.qinvmR[0]);
  perf_intt(p, r_1, rns_params.q[0], rns_params.R_q[0], rns_params.qinvmR[0]);

  // Perf CKKS Operations
  perf_keygen(&ckks_params);
  ciphertext_new(&ct, &ckks_params.rns);
  perf_encrypt(&ct, p, &v, &ckks_params);
  perf_decrypt(&v, p, &ct, &ckks_params);

  // Free
  vec_free(&v);
  vec_free(&w);
  poly_free(&p);
  poly_free(&r);
  poly_free(&r_1);
  rns_params_free(&rns_params);
  ckks_free(&ckks_params);
  ciphertext_free(&ct);
}
