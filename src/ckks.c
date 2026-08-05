
#include "ckks.h"

#include <string.h>

void ckks_init(ckks_params* params) {
  gen_w(params->w, params->n);
  poly s_hat;
  poly_new(&s_hat, params->n);
  poly_ternary_sampler(s_hat, params->n);
  for (uint32_t i = 0; i < params->l; i++) {
    gen_eta(params->eta[i], params->n, params->q[i]);
    memcpy(params->s[i], s_hat, sizeof(int64_t) * params->n);
    poly_ntt(params->s[i], params->eta[i], params->n);
  }
  uint8_t seed[SEED_BYTES];
  randombytes(seed, SEED_BYTES);
  // Consider that the p values are already on the NTT domain
  for (uint32_t i = 0; i < params->l; i++) {
    poly_uniform_sampler(params->p[i], seed, params->q[i], params->n);
    seed[0] += 1;
  }
}

void ckks_encode(poly m, const vec* v, const vec* w, int64_t scale,
                 uint32_t n) {
  // First step: Convert v into hermitian v2 = (v, v)
  vec v2;
  vec_new(&v2, n, v->int_part);
  vec_hermitian(&v2, v, n);
  // Second step: IFFT
  ifft(&v2, w, n);
  // Third step: Rescale
  vector_to_poly(m, &v2, scale, n);
  // Free
  vec_free(&v2);
}

void ckks_decode(vec* v, const poly m, const vec* w, int64_t scale,
                 uint32_t n) {
  // First step: Rescale
  vec v2;
  vec_new(&v2, n, v->int_part);
  poly_to_vector(&v2, m, scale, n);
  // Second step: FFT
  fft(&v2, w, n);
  // Third step: Convert into input complex vector
  vec_input(v, &v2, n);
  // Free
  vec_free(&v2);
}

void ckks_encrypt(rlwe ct, const poly m, const poly s, const int64_t* w,
                  uint32_t n, const int64_t* q, const uint32_t l) {
  poly e;
  poly_new(&e, n);
  // DeltaM + E
  poly_cbd24_sampler(e, n);
  poly_add(e, m, n);
  // -AS + DeltaM + E
  for (uint32_t i = 0; i < l; i++) {
    poly_uniform_sampler(ct.b[i], q[i], n);
    poly_ntt(s[i], w, n);
  }

  poly_free(&e);
}
