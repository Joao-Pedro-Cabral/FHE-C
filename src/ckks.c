
#include "ckks.h"

#include <stdlib.h>
#include <string.h>

#include "mod.h"
#include "random.h"
#include "utils.h"

void ckks_keygen(ckks_params_t* params, uint32_t q0_bits, uint32_t scale_bits,
                 uint32_t n2, uint32_t l, uint32_t n, int set_n,
                 uint64_t* q_chain, int set_q) {
  // RNS construction
  rns_params_new(&params->rns, q0_bits, scale_bits, l, n, set_n, q_chain,
                 set_q);
  // complex roots
  params->n2 = n2;
  vec_new(&params->w, n2);
  // XOF state
  uint8_t seed[SEED_BYTES];
  randombytes(seed, SEED_BYTES);
  xof_absorb(&params->xof_state, seed);
  // keys
  rlwe_keys_new(&params->keys, &params->rns);
#pragma omp parallel
  {
    vec_gen_roots(&params->w, n2);
    rlwe_gen_secret_key(&params->keys, &params->xof_state, &params->rns);
  }
}

void ckks_relinearization_keygen(ckks_params_t* params) {
  rlwe_gen_relinearization_key(&params->keys, &params->xof_state, &params->rns);
}

// KeySwitch(S^{-1})
void ckks_conjugation_keygen(ckks_params_t* params) {
  poly_t s2;
  poly_new(&s2, params->rns.n);
#pragma omp parallel
  {
    poly_t e, me;
    poly_new(&e, params->rns.n);
    poly_new(&me, params->rns.n);
    uint32_t exp = 2 * params->rns.n - 1;
    poly_permute_omp(s2, params->keys.s_hat, exp, params->rns.n);
    rlwe_key_switch_encrypt_omp(params->keys.cg, &params->xof_state, s2, me, e,
                                &params->keys, &params->rns);
    poly_free(&e);
    poly_free(&me);
  }
  poly_free(&s2);
}

// RLWE(0)
void ckks_public_keygen(ckks_params_t* params) {
  poly_t m;
  poly_new(&m, params->rns.n);
  memset(m, 0, sizeof(int64_t) * params->rns.n);
  rlwe_plaintext_encrypt(&params->keys.pk, &params->xof_state, m, &params->keys,
                         &params->rns);
  poly_free(&m);
}

// Trade-off memory/performance for inner-product:
// (log2_n - 1) keys + (log2_n - 1) rotations for reduction
// log2_n - 1 because we work with n/2 complex vector slots
// Our rotation keys are 5^(n/4), 5^(n/8), ... 5^1
// We use these keys to reduce the memory usage during reductions
void ckks_inner_product_keygen(ckks_params_t* params) {
  uint32_t log2_n2 = floor_log2(params->rns.n) - 1;
  uint64_t max_n = 0;
#pragma omp parallel
  {
    poly_t s2, e, me;
    poly_new(&s2, params->rns.n);
    poly_new(&e, params->rns.n);
    poly_new(&me, params->rns.n);
#pragma omp for reduction(max : max_n)
    for (uint32_t i = 0; i < log2_n2; i++) {
      xof_state_t xof_state2;
      memcpy(&xof_state2, &params->xof_state, sizeof(xof_state2));
      uint32_t exp = params->rns.n >> (i + 2);
      uint32_t g = mod_pow2(5, exp, 2 * params->rns.n);
      poly_permute(s2, params->keys.s_hat, g, params->rns.n);
      rlwe_key_switch_encrypt(params->keys.ip[i], &xof_state2, s2, me, e,
                              &params->keys, &params->rns);
      max_n = (xof_state2.n > max_n) ? xof_state2.n : max_n;
    }
    poly_free(&s2);
    poly_free(&e);
    poly_free(&me);
  }
  xof_state_update(&params->xof_state, max_n);
}

// For linear transformations, we need n1+n2 keys for Dual Hoisting:
// n1: 0, 1, ..., n1-1
// n2: 0*n1, 1*n1, ..., (n2-1)*n1
// Where n1*n2 = n/2 = m and n1 = N1N2_RATIO*n2
void ckks_linear_transform_keygen(ckks_params_t* params) {
  uint32_t log2_m = floor_log2(params->rns.n) - 1;
  uint32_t log2_n2 = floor_log2(params->rns.n / (2 * N1N2_RATIO)) / 2;
  uint32_t log2_n1 = log2_m - log2_n2;
  uint32_t n1 = (1 << log2_n1);
  uint32_t n2 = (1 << log2_n2);
  uint64_t max_n = 0;
#pragma omp parallel
  {
    poly_t s2, e, me;
    poly_new(&s2, params->rns.n);
    poly_new(&e, params->rns.n);
    poly_new(&me, params->rns.n);
#pragma omp for reduction(max : max_n)
    for (uint32_t i = 0; i < n1 + n2; i++) {
      xof_state_t xof_state2;
      memcpy(&xof_state2, &params->xof_state, sizeof(xof_state2));
      uint32_t exp = (i >= n1) ? n1 * (i - n1) : i;
      uint32_t g = mod_pow2(5, exp, 2 * params->rns.n);
      poly_permute(s2, params->keys.s_hat, g, params->rns.n);
      rlwe_key_switch_encrypt(params->keys.lt[i], &xof_state2, s2, me, e,
                              &params->keys, &params->rns);
      max_n = (xof_state2.n > max_n) ? xof_state2.n : max_n;
    }
    poly_free(&s2);
    poly_free(&e);
    poly_free(&me);
  }
  xof_state_update(&params->xof_state, max_n);
}

void ckks_free(ckks_params_t* params) {
  vec_free(&params->w);
  rns_params_free(&params->rns);
  rlwe_keys_free(&params->keys, &params->rns);
}

void ckks_encode(poly_t m, const vec_t* v, const ckks_params_t* params) {
  vec_t v2;
  vec_new(&v2, params->n2);
  memcpy(v2.rel, v->rel, 2 * params->n2 * sizeof(double));
#pragma omp parallel
  {
    // First step: IFFT
    ifft(&v2, &params->w, params->n2);
    // Second step: Rescale (with base scale) + Bit Reverse
    vector_to_poly(m, &v2, params->rns.scales[params->rns.l - 1], params->n2,
                   params->rns.n);
  }
  // Free
  vec_free(&v2);
}

void ckks_decode(vec_t* v, const poly_t m, uint32_t level,
                 const ckks_params_t* params) {
#pragma omp parallel
  {
    // First step: Rescale (accordingly to the current level) + Bit Reverse
    poly_to_vector(v, m, params->rns.scales[level - 1], params->n2,
                   params->rns.n);
    // Second step: FFT
    fft(v, &params->w, params->n2);
  }
}

void ckks_encrypt(ciphertext_t* ct, const poly_t m, ckks_params_t* params) {
  rlwe_plaintext_encrypt(&ct->ab, &params->xof_state, m, &params->keys,
                         &params->rns);
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

void ckks_server_encrypt(ciphertext_t* ct, const poly_t m,
                         ckks_params_t* params) {
  rlwe_plaintext_encrypt_server(&ct->ab, &params->xof_state, m, params->keys.pk,
                                &params->rns);
}

void ckks_ciphertext_conjugation(ciphertext_t* ct2, const ciphertext_t* ct1,
                                 const ckks_params_t* params) {
  rns_t tmp_r, tmp_r2;
  poly_t tmp_p, tmp_p2;
  uint32_t n = params->rns.n;
  rns_new(&tmp_r, &params->rns);
  rns_new(&tmp_r2, &params->rns);
  poly_new(&tmp_p, n);
  poly_new(&tmp_p2, n);
  ciphertext_rotation(ct2, ct1, 2 * params->rns.n - 1, tmp_r, tmp_r2, tmp_p,
                      tmp_p2, params->keys.cg, &params->rns);
  rns_free(&tmp_r);
  rns_free(&tmp_r2);
  poly_free(&tmp_p);
  poly_free(&tmp_p2);
}

void ckks_ciphertext_inner_product(ciphertext_t* ct3, const ciphertext_t* ct1,
                                   const ciphertext_t* ct2,
                                   const ckks_params_t* params) {
  ciphertext_t tmp_c;
  rns_t tmp_r, tmp_r2;
  poly_t tmp_p, tmp_p2;
  uint32_t n = params->rns.n;
  ciphertext_new(&tmp_c, &params->rns);
  rns_new(&tmp_r, &params->rns);
  rns_new(&tmp_r2, &params->rns);
  poly_new(&tmp_p, n);
  poly_new(&tmp_p2, n);
  ciphertext_ciphertext_multiplication(ct3, ct1, ct2, params->keys.r,
                                       &params->rns);
  uint32_t log2_n2 = floor_log2(n) - 1;
  for (uint32_t i = 0; i < log2_n2; i++) {
    uint32_t exp = n >> (i + 2);
    exp = mod_pow2(5, exp, 2 * n);
    ciphertext_rotation(&tmp_c, ct3, exp, tmp_r, tmp_r2, tmp_p, tmp_p2,
                        params->keys.ip[i], &params->rns);
    ciphertext_ciphertext_addition(ct3, ct3, &tmp_c, &params->rns);
  }
  ciphertext_free(&tmp_c);
  rns_free(&tmp_r);
  rns_free(&tmp_r2);
  poly_free(&tmp_p);
  poly_free(&tmp_p2);
}

void ckks_ciphertext_linear_transform(ciphertext_t* ct2,
                                      const ciphertext_t* ct1,
                                      const poly_t* mat,
                                      const ckks_params_t* params) {
  ciphertext_t acc, tmp_c;
  rns_t tmp_r, tmp_r2;
  poly_t tmp_p, tmp_p2;
  uint32_t n = params->rns.n;
  uint32_t log2_m = floor_log2(n) - 1;
  uint32_t log2_n2 = floor_log2(n / (2 * N1N2_RATIO)) / 2;
  uint32_t log2_n1 = log2_m - log2_n2;
  uint32_t n1 = (1 << log2_n1);
  uint32_t n2 = (1 << log2_n2);
  ciphertext_new(&acc, &params->rns);
  ciphertext_new(&tmp_c, &params->rns);
  rns_new(&tmp_r, &params->rns);
  rns_new(&tmp_r2, &params->rns);
  poly_new(&tmp_p, n);
  poly_new(&tmp_p2, n);
  // Compute rot(ct1, i); i = 0, 1, ..., n1-1 (jump 0)
  ciphertext_t* rot_ct1 = malloc(n1 * sizeof(ciphertext_t));
  for (uint32_t i = 1; i < n1; i++) {
    ciphertext_new(&rot_ct1[i], &params->rns);
    uint32_t exp = mod_pow2(5, i, 2 * n);
    ciphertext_rotation(&rot_ct1[i], ct1, exp, tmp_r, tmp_r2, tmp_p, tmp_p2,
                        params->keys.lt[i], &params->rns);
  }
  // Zero final ciphertext
  ct2->level = ct1->level - 1;
  for (uint32_t i = 0; i < ct2->level; i++) {
    memset(&ct2->ab.a[n * i], 0, n * sizeof(int64_t));
    memset(&ct2->ab.b[n * i], 0, n * sizeof(int64_t));
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
    ciphertext_mod_switch(&acc, tmp_r, &params->rns);
    // Rotate and Accumulate
    uint32_t exp = mod_pow2(5, n1 * i, 2 * n);
    ciphertext_rotation(&tmp_c, &acc, exp, tmp_r, tmp_r2, tmp_p, tmp_p2,
                        params->keys.lt[n1 + i], &params->rns);
    ciphertext_ciphertext_addition(ct2, ct2, &tmp_c, &params->rns);
  }

  for (uint32_t i = 1; i < n1; i++) {
    ciphertext_free(&rot_ct1[i]);
  }
  free(rot_ct1);
  ciphertext_free(&tmp_c);
  ciphertext_free(&acc);
  rns_free(&tmp_r);
  rns_free(&tmp_r2);
  poly_free(&tmp_p);
  poly_free(&tmp_p2);
}

// Encode a linear transformation with the correct rotations
void ckks_linear_transform_encode(poly_t* mat2, const vec_t* mat1,
                                  const ckks_params_t* params) {
  vec_t diag;
  uint32_t m = params->n2;
  uint32_t n = params->rns.n;
  uint32_t ratio = n / (2 * m);
  // Zero diagonals for sparse encoding
  for (uint32_t i = 0; i < n; i++) {
    if (i % ratio != 0) {
      memset(mat2[i], 0, sizeof(int64_t) * n);
    }
  }
  uint32_t log2_m = floor_log2(m);
  uint32_t log2_n2 = floor_log2(m / (N1N2_RATIO)) / 2;
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
    vector_to_poly(mat2[i * ratio], &diag,
                   params->rns.scales[params->rns.l - 1], m, params->rns.n);
  }
  vec_free(&diag);
}
