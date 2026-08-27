
#include "ckks.h"

#include <stdlib.h>
#include <string.h>

#include "random.h"
#include "utils.h"

void ckks_keygen(ckks_params_t* params, uint32_t q0_bits, uint32_t scale_bits,
                 uint32_t l) {
  // RNS construction
  rns_params_new(&params->rns, q0_bits, scale_bits, l);
  // complex roots
  vec_new(&params->w, params->rns.n / 2);
  vec_gen_roots(&params->w, params->rns.n / 2);
  // secret key for RLWE
  rlwe_keys_new(&params->keys, &params->rns);
  rlwe_gen_secret_key(&params->keys, &params->rns);
}

void ckks_relinearization_keygen(ckks_params_t* params) {
  rlwe_gen_relinearization_key(&params->keys, &params->rns);
}

// Trade-off memory/performance for inner-product:
// (log2_n - 1) keys + (log2_n - 1) rotations for reduction
// log2_n - 1 because we work with n/2 complex vector slots
// Our rotation keys are 5^(n/4), 5^(n/8), ... 5^1
// We use these keys to reduce the memory usage during reductions
void ckks_inner_product_keygen(ckks_params_t* params) {
  poly_t s2, e, me;
  poly_new(&s2, params->rns.n);
  poly_new(&e, params->rns.n);
  poly_new(&me, params->rns.n);
  uint32_t log2_n2 = floor_log2(params->rns.n) - 1;
  for (uint32_t i = 0; i < log2_n2; i++) {
    uint32_t exp = params->rns.n >> (i + 2);
    uint32_t g = mod_pow2(5, exp, 2 * params->rns.n);
    poly_permute(s2, params->keys.s_hat, g, params->rns.n);
    rlwe_key_switch_encrypt(params->keys.ip[i], s2, me, e, &params->keys,
                            &params->rns);
  }
  poly_free(&s2);
  poly_free(&e);
  poly_free(&me);
}

// For linear transformations, we need n1+n2 keys for Dual Hoisting:
// n1: 0, 1, ..., n1-1
// n2: 0*n1, 1*n1, ..., (n2-1)*n1
// Where n1*n2 = n/2 = m and n1 = N1N2_RATIO*n2
void ckks_linear_transform_keygen(ckks_params_t* params) {
  poly_t s2, e, me;
  poly_new(&s2, params->rns.n);
  poly_new(&e, params->rns.n);
  poly_new(&me, params->rns.n);
  uint32_t log2_m = floor_log2(params->rns.n) - 1;
  uint32_t log2_n2 = floor_log2(params->rns.n / (2 * N1N2_RATIO)) / 2;
  uint32_t log2_n1 = log2_m - log2_n2;
  uint32_t n1 = (1 << log2_n1);
  uint32_t n2 = (1 << log2_n2);
  for (uint32_t i = 0; i < n1; i++) {
    uint32_t exp = i;
    uint32_t g = mod_pow2(5, exp, 2 * params->rns.n);
    poly_permute(s2, params->keys.s_hat, g, params->rns.n);
    rlwe_key_switch_encrypt(params->keys.lt[i], s2, me, e, &params->keys,
                            &params->rns);
  }
  for (uint32_t i = 0; i < n2; i++) {
    uint32_t exp = n1 * i;
    uint32_t g = mod_pow2(5, exp, 2 * params->rns.n);
    poly_permute(s2, params->keys.s_hat, g, params->rns.n);
    rlwe_key_switch_encrypt(params->keys.lt[n1 + i], s2, me, e, &params->keys,
                            &params->rns);
  }
  poly_free(&s2);
  poly_free(&e);
  poly_free(&me);
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
  // TODO: PLEASE DO A FASTER INNER PRODUCT WITHOUT CT4
  ciphertext_t tmp_c, ct4;
  rns_t tmp_r;
  poly_t tmp_p;
  ciphertext_new(&tmp_c, &params->rns);
  ciphertext_new(&ct4, &params->rns);
  rns_new(&tmp_r, &params->rns);
  poly_new(&tmp_p, params->rns.n);
  ciphertext_ciphertext_multiplication(ct3, ct1, ct2, params->keys.r,
                                       &params->rns);
  ct4.level = ct3->level;
  uint32_t log2_n2 = floor_log2(params->rns.n) - 1;
  for (uint32_t i = 0; i < log2_n2; i++) {
    for (uint32_t k = 0; k < ct4.level; k++) {
      memcpy(ct4.ab.a[k], ct3->ab.a[k], params->rns.n * sizeof(int64_t));
      memcpy(ct4.ab.b[k], ct3->ab.b[k], params->rns.n * sizeof(int64_t));
    }
    uint32_t exp = params->rns.n >> (i + 2);
    ciphertext_rotation(&tmp_c, &ct4, exp, tmp_r, tmp_p, params->keys.ip[i],
                        &params->rns);
    ciphertext_ciphertext_addition(ct3, ct3, &tmp_c, &params->rns);
  }
  ciphertext_free(&tmp_c, &params->rns);
  ciphertext_free(&ct4, &params->rns);
  rns_free(&tmp_r, &params->rns);
  poly_free(&tmp_p);
}

void ckks_ciphertext_linear_transform(ciphertext_t* ct2,
                                      const ciphertext_t* ct1,
                                      const poly_t* mat,
                                      const ckks_params_t* params) {
  ciphertext_t acc, tmp_c;
  rns_t tmp_r;
  poly_t tmp_p;
  uint32_t log2_m = floor_log2(params->rns.n) - 1;
  uint32_t log2_n2 = floor_log2(params->rns.n / (2 * N1N2_RATIO)) / 2;
  uint32_t log2_n1 = log2_m - log2_n2;
  uint32_t n1 = (1 << log2_n1);
  uint32_t n2 = (1 << log2_n2);
  ciphertext_new(&acc, &params->rns);
  ciphertext_new(&tmp_c, &params->rns);
  rns_new(&tmp_r, &params->rns);
  poly_new(&tmp_p, params->rns.n);
  // Compute rot(ct1, i); i = 0, 1, ..., n1-1 (jump 0)
  ciphertext_t* rot_ct1 = malloc(n1 * sizeof(ciphertext_t));
  for (uint32_t i = 1; i < n1; i++) {
    ciphertext_new(&rot_ct1[i], &params->rns);
    ciphertext_rotation(&rot_ct1[i], ct1, i, tmp_r, tmp_p, params->keys.lt[i],
                        &params->rns);
  }
  // Zero final ciphertext
  ct2->level = ct1->level - 1;
  for (uint32_t i = 0; i < ct2->level; i++) {
    memset(ct2->ab.a[i], 0, params->rns.n * sizeof(int64_t));
    memset(ct2->ab.b[i], 0, params->rns.n * sizeof(int64_t));
  }
  for (uint32_t i = 0; i < n2; i++) {
    ciphertext_plaintext_base_multiplication(&acc, ct1, mat[n1 * i],
                                             &params->rns);
    // Multiply and accumulate without mod switch
    for (uint32_t j = 1; j < n1; j++) {
      ciphertext_plaintext_base_multiplication(&tmp_c, &rot_ct1[j],
                                               mat[n1 * i + j], &params->rns);
      ciphertext_ciphertext_addition(&acc, &acc, &tmp_c, &params->rns);
    }
    ciphertext_mod_switch(&acc, tmp_p, &params->rns);
    // Rotate and Accumulate
    ciphertext_rotation(&tmp_c, &acc, n1 * i, tmp_r, tmp_p,
                        params->keys.lt[n1 + i], &params->rns);
    ciphertext_ciphertext_addition(ct2, ct2, &tmp_c, &params->rns);
  }

  for (uint32_t i = 1; i < n1; i++) {
    ciphertext_free(&rot_ct1[i], &params->rns);
  }
  free(rot_ct1);
  ciphertext_free(&tmp_c, &params->rns);
  ciphertext_free(&acc, &params->rns);
  rns_free(&tmp_r, &params->rns);
  poly_free(&tmp_p);
}

// Encode a linear transformation with the correct rotations
void ckks_linear_transform_encode(poly_t* mat2, const vec_t* mat1,
                                  const ckks_params_t* params) {
  vec_t diag;
  uint32_t m = params->rns.n / 2;
  uint32_t log2_m = floor_log2(m);
  uint32_t log2_n2 = floor_log2(params->rns.n / (2 * N1N2_RATIO)) / 2;
  uint32_t log2_n1 = log2_m - log2_n2;
  uint32_t n1 = (1 << log2_n1);
  vec_new(&diag, m);
  for (uint32_t i = 0; i < m; i++) {
    uint32_t n1_i = (i / n1) * n1;
    // Create rotated diagonal by n1*i to the right
    for (uint32_t j = 0; j < m; j++) {
      diag.rel[(n1_i + j) & (m - 1)] = mat1[j].rel[(i + j) & (m - 1)];
      diag.img[(n1_i + j) & (m - 1)] = mat1[j].img[(i + j) & (m - 1)];
    }
    // Encode
    ifft(&diag, &params->w, m);
    vector_to_poly(mat2[i], &diag, params->rns.scales[params->rns.l - 1],
                   params->rns.n);
  }
  vec_free(&diag);
}
