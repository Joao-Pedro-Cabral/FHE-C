
#include "ckks.h"

#include <string.h>

#include "random.h"
#include "utils.h"

void ckks_keygen(ckks_params_t* params, uint32_t q0_bits, uint32_t scale_bits,
                 uint32_t l) {
  // RNS construction
  rns_params_new(&params->rns, q0_bits, scale_bits, l);
  // complex roots
  vec_new(&params->w, params->rns.n);
  vec_gen_roots(&params->w, params->rns.n);
  // secret key for RLWE
  rlwe_keys_new(&params->keys, &params->rns);
  rlwe_gen_secret_key(&params->keys, &params->rns);
}

void ckks_relinearization_keygen(ckks_params_t* params) {
  rlwe_gen_relinearization_key(&params->keys, &params->rns);
}

void ckks_rotation_keygen(ckks_params_t* params) {
  rlwe_gen_rotation_key(&params->keys, &params->rns);
}

void ckks_free(ckks_params_t* params) {
  vec_free(&params->w);
  rns_params_free(&params->rns);
  rlwe_keys_free(&params->keys, &params->rns);
}

void ckks_encode(poly_t m, const vec_t* v, const ckks_params_t* params) {
  vec_t v2;
  vec_new(&v2, params->rns.n / 2);
  memcpy(v2.rel, v->rel, params->rns.n * sizeof(double));
  // First step: IFFT
  ifft(&v2, &params->w, params->rns.n / 2);
  // Second step: Rescale (with base scale) + Bit Reverse
  vector_to_poly(m, &v2, params->rns.scales[params->rns.l - 1], params->rns.n);
  // Free
  vec_free(&v2);
}

void ckks_decode(vec_t* v, const poly_t m, uint32_t level,
                 const ckks_params_t* params) {
  // First step: Rescale (accordingly to the current level) + Bit Reverse
  poly_to_vector(v, m, params->rns.scales[level - 1], params->rns.n);
  // Second step: FFT
  fft(v, &params->w, params->rns.n / 2);
}

void ckks_encrypt(ciphertext_t* ct, const poly_t m,
                  const ckks_params_t* params) {
  rlwe_plaintext_encrypt(&ct->ab, m, &params->keys, &params->rns);
}

void ckks_decrypt(poly_t m, const ciphertext_t* ct,
                  const ckks_params_t* params) {
  rlwe_plaintext_decrypt(m, &ct->ab, &params->keys, &params->rns);
}

void ckks_ciphertext_ciphertext_addition(ciphertext_t* ct3,
                                         const ciphertext_t* ct1,
                                         const ciphertext_t* ct2,
                                         const ckks_params_t* params) {
  ciphertext_ciphertext_addition(ct3, ct1, ct2, &params->rns);
}

void ckks_ciphertext_ciphertext_multiplication(ciphertext_t* ct3,
                                               const ciphertext_t* ct1,
                                               const ciphertext_t* ct2,
                                               const ckks_params_t* params) {
  ciphertext_ciphertext_multiplication(ct3, ct1, ct2, params->keys.r,
                                       &params->rns);
}

void ckks_ciphertext_plaintext_addition(ciphertext_t* ct2,
                                        const ciphertext_t* ct1, const poly_t m,
                                        const ckks_params_t* params) {
  ciphertext_plaintext_addition(ct2, ct1, m, &params->rns);
}

void ckks_ciphertext_plaintext_multiplication(ciphertext_t* ct2,
                                              const ciphertext_t* ct1,
                                              const poly_t m,
                                              const ckks_params_t* params) {
  ciphertext_plaintext_multiplication(ct2, ct1, m, &params->rns);
}

void ckks_ciphertext_inner_product(ciphertext_t* ct3, const ciphertext_t* ct1,
                                   const ciphertext_t* ct2,
                                   const ckks_params_t* params) {
  ciphertext_t tmp;
  rns_t a;
  poly_t tmp2;
  ciphertext_new(&tmp, &params->rns);
  rns_new(&a, &params->rns);
  poly_new(&tmp2, params->rns.n);
  ciphertext_ciphertext_multiplication(ct3, ct1, ct2, params->keys.r,
                                       &params->rns);
  uint32_t log2_n2 = floor_log2(params->rns.n) - 1;
  uint32_t exp = 0;
  for (uint32_t i = 0; i < log2_n2; i++) {
    exp += params->rns.n >> (i + 2);
    ciphertext_rotation(&tmp, ct3, exp, a, tmp2, params->keys.h[i],
                        &params->rns);
    ciphertext_ciphertext_addition(ct3, ct3, &tmp, &params->rns);
  }
  ciphertext_free(&tmp, &params->rns);
  rns_free(&a, &params->rns);
  poly_free(&tmp2);
}
