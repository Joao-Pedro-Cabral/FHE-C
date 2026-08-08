
#include "ckks.h"

#include "random.h"

void ckks_new(ckks_params* params, uint32_t n, uint32_t l, const uint64_t* q,
              int64_t scale) {
  params->n = n;
  params->l = l;
  params->q = q;
  params->scale = scale;
  // w
  vec_new(&params->w, n);
  vec_gen_w(&params->w, n);
  // eta for NTT
  set eta0;
  set_new(&eta0, n, 1);
  set_new(&params->eta, n, l);
  set_gen_eta(params->eta, eta0[0], q, n, l);
  // s
  rlwe_new(&params->s, n, l);
  rlwe_gen_secret(params->s, params->eta, q, n, l);
  // eta for INTT
  set_gen_reverse_eta(params->eta, eta0[0], q, n, l);
  set_free(&eta0, 1);
}

void ckks_free(ckks_params* params) {
  vec_free(&params->w);
  set_free(&params->eta, params->l);
  rlwe_free(&params->s, params->l);
}

void ckks_encode(poly m, const vec* v, const ckks_params* params) {
  // First step: Convert v into hermitian v2 = (v, v)
  vec v2;
  vec_new(&v2, params->n);
  vec_hermitian(&v2, v, params->n);
  // Second step: IFFT
  ifft(&v2, &params->w, params->n);
  // Third step: Rescale
  vector_to_poly(m, &v2, params->scale, params->n);
  // Free
  vec_free(&v2);
}

void ckks_decode(vec* v, const poly m, const ckks_params* params) {
  // First step: Rescale
  vec v2;
  vec_new(&v2, params->n);
  poly_to_vector(&v2, m, params->scale, params->n);
  // Second step: FFT
  fft(&v2, &params->w, params->n);
  // Third step: Convert into input complex vector
  vec_input(v, &v2, params->n);
  // Free
  vec_free(&v2);
}

void ckks_encrypt(cyphertext* ct, const poly m, const ckks_params* params) {
  // DeltaM + E
  poly e;
  poly_new(&e, params->n);
  poly_cbd24_sampler(e, params->n);
  poly_add_no_mod(e, m, params->n);
  // Init A
  randombytes(ct->seed, SEED_BYTES);
  // -AS + DeltaM + E
  for (uint32_t i = 0; i < params->l; i++) {
    poly_uniform_sampler(ct->b[i], ct->seed, params->q[i], params->n);
    ct->seed[0] += 1;
    poly_pointwise_mul(ct->b[i], params->s[i], params->q[i], params->n);
    poly_intt(ct->b[i], params->eta[i], params->q[i], params->n);
    poly_rsub(ct->b[i], e, params->q[i], params->n);
  }
  ct->seed[0] -= (uint8_t)params->l;  // TODO: l < 255
  poly_free(&e);
}

void ckks_decrypt(poly m, const cyphertext* ct, const ckks_params* params) {
  // For decrypting we only need l = 0 (deltaM + E < q_0/2)
  // B + AS = DeltaM + E
  poly_uniform_sampler(m, ct->seed, params->q[0], params->n);
  poly_pointwise_mul(m, params->s[0], params->q[0], params->n);
  poly_intt(m, params->eta[0], params->q[0], params->n);
  poly_add(m, ct->b[0], params->q[0], params->n);
}

void cyphertext_new(cyphertext* ct, const ckks_params* params) {
  rlwe_new(&(ct->b), params->n, params->l);
}

void cyphertext_free(cyphertext* ct, const ckks_params* params) {
  rlwe_free(&(ct->b), params->l);
}
