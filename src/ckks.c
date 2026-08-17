
#include "ckks.h"

#include "random.h"

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
  // First step: Convert v into hermitian v2 = (v, v)
  vec_t v2;
  vec_new(&v2, params->rns.n);
  vec_hermitian(&v2, v, params->rns.n);
  // Second step: IFFT
  ifft(&v2, &params->w, params->rns.n);
  // Third step: Rescale (with base scale)
  vector_to_poly(m, &v2, params->rns.scales[params->rns.l - 1], params->rns.n);
  // Free
  vec_free(&v2);
}

void ckks_decode(vec_t* v, const poly_t m, uint32_t level,
                 const ckks_params_t* params) {
  // First step: Rescale (accordingly to the current level)
  vec_t v2;
  vec_new(&v2, params->rns.n);
  poly_to_vector(&v2, m, params->rns.scales[level - 1], params->rns.n);
  // Second step: FFT
  fft(&v2, &params->w, params->rns.n);
  // Third step: Convert into input complex vector
  vec_input(v, &v2, params->rns.n);
  // Free
  vec_free(&v2);
}

void ckks_encrypt(ciphertext_t* ct, const poly_t m,
                  const ckks_params_t* params) {
  rlwe_plaintext_encrypt(&ct->ab, m, &params->keys, &params->rns);
}

void ckks_decrypt(poly_t m, const ciphertext_t* ct,
                  const ckks_params_t* params) {
  rlwe_plaintext_decrypt(m, &ct->ab, &params->keys, &params->rns);
}

void ckks_ciphertext_ciphertext_multiplication(ciphertext_t* ct3,
                                               const ciphertext_t* ct1,
                                               const ciphertext_t* ct2,
                                               const ckks_params_t* params) {
  ciphertext_ciphertext_multiplication(ct3, ct1, ct2, params->keys.r,
                                       &params->rns);
}
