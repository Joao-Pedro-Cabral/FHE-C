
#include "ckks.h"

#include <string.h>
#include "random.h"

void ckks_init(ckks_params* params) {
  gen_w(params->w, params->n);
  poly s_hat;
  poly_new(&s_hat, params->n);
  poly_ternary_sampler(s_hat, params->n);
  // for (uint32_t i = 0; i < params->l; i++) {
  //   gen_eta(params->eta[i], params->n, params->q[i]);
  //   memcpy(params->s[i], s_hat, sizeof(int64_t) * params->n);
  //   poly_ntt(params->s[i], params->eta[i], params->n);
  // }
  poly_free(&s_hat);
}

void ckks_encode(poly m, const vec* v, const ckks_params* params) {
  // First step: Convert v into hermitian v2 = (v, v)
  vec v2;
  vec_new(&v2, params->n, v->int_part);
  vec_hermitian(&v2, v, params->n);
  // Second step: IFFT
  ifft(&v2, params->w, params->n);
  // Third step: Rescale
  vector_to_poly(m, &v2, params->scale, params->n);
  // Free
  vec_free(&v2);
}

void ckks_decode(vec* v, const poly m, const ckks_params* params) {
  // First step: Rescale
  vec v2;
  vec_new(&v2, params->n, v->int_part);
  poly_to_vector(&v2, m, params->scale, params->n);
  // Second step: FFT
  fft(&v2, params->w, params->n);
  // Third step: Convert into input complex vector
  vec_input(v, &v2, params->n);
  // Free
  vec_free(&v2);
}

void ckks_encrypt(cyphertext ct, const poly m, const ckks_params* params) {
  // DeltaM + E
  poly e;
  poly_new(&e, params->n);
  poly_cbd24_sampler(e, params->n);
  poly_add(e, m, params->n);
  // Init A
  randombytes(ct.seed, SEED_BYTES);
  // -AS + DeltaM + E
  for (uint32_t i = 0; i < params->l; i++) {
    poly_uniform_sampler(ct.b[i], ct.seed, params->q[i], params->n);
    ct.seed[0] += 1;
    poly_pointwise_mul(ct.b[i], params->s[i], params->n);
    poly_rsub(ct.b[i], e, params->n);
  }
  ct.seed[0] -= (uint8_t)params->l; // TODO: l < 255
  poly_free(&e);
}

void ckks_decrypt(poly m, const cyphertext ct, const ckks_params* params) {
  // For decrypting we only need l = 0 (deltaM + E < q_0/2)
  // B + AS = DeltaM + E
  poly_uniform_sampler(m, ct.seed, params->q[0], params->n);
  poly_pointwise_mul(m, params->s[0], params->n);
  poly_add(m, ct.b[0], params->n);
}
