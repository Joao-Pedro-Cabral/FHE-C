
#include "ckks.h"

#include "random.h"

void ckks_keygen(ckks_params_t* params, uint32_t q0_bits, uint32_t scale_bits,
                 uint32_t l, uint64_t scale) {
  params->scale = scale;
  // RNS construction
  rns_params_new(&params->rns, q0_bits, scale_bits, l);
  // complex roots
  vec_new(&params->w, params->rns.n);
  vec_gen_roots(&params->w, params->rns.n);
  // roots for RLWE
  rns_gen_roots(&params->rns);
  // secret key for RLWE
  rlwe_keys_new(&params->keys, &params->rns);
  rlwe_gen_secret_key(&params->keys, &params->rns);
  // roots for INTT
  rns_gen_inverse_roots(&params->rns);
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
  // Third step: Rescale
  vector_to_poly(m, &v2, params->scale, params->rns.n);
  // Free
  vec_free(&v2);
}

void ckks_decode(vec_t* v, const poly_t m, const ckks_params_t* params) {
  // First step: Rescale
  vec_t v2;
  vec_new(&v2, params->rns.n);
  poly_to_vector(&v2, m, params->scale, params->rns.n);
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

void cyphertext_new(ciphertext_t* ct, const ckks_params_t* params) {
  rlwe_new(&ct->ab, &params->rns);
  ct->n = params->rns.n;
  ct->l = params->rns.l;
  ct->scale = params->scale;
  ct->int_scale = 1;
  ct->noise_scale = 1;
  ct->encoding_type = 4;  // OpenFHE: CKKS_PACKED_ENCODING
}

void cyphertext_free(ciphertext_t* ct, const ckks_params_t* params) {
  rlwe_free(&ct->ab, &params->rns);
}
