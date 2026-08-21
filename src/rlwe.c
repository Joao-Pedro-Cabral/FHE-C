
#include "rlwe.h"

#include <stdlib.h>
#include <string.h>

#include "mod.h"
#include "utils.h"

void rlwe_new(rlwe_t* ab, const rns_params_t* params) {
  rns_new(&ab->b, params);
  rns_new(&ab->a, params);
}

void rlwe_free(rlwe_t* ab, const rns_params_t* params) {
  rns_free(&ab->b, params);
  rns_free(&ab->a, params);
}

void rlwe_keys_new(keys_t* keys, const rns_params_t* params) {
  poly_new(&keys->s_hat, params->n);
  rns_new(&keys->s, params);
  keys->r = malloc(params->l * sizeof(rlwe_t));
  for (uint32_t i = 0; i < params->l; i++) {
    rlwe_new(&keys->r[i], params);
  }
  uint32_t log2_n = floor_log2(params->n) - 1;
  keys->h = malloc(log2_n * sizeof(switch_t));
  for (uint32_t i = 0; i < log2_n; i++) {
    keys->h[i] = malloc(params->l * sizeof(rlwe_t));
    for (uint32_t j = 0; j < params->l; j++) {
      rlwe_new(&keys->h[i][j], params);
    }
  }
}

void rlwe_keys_free(keys_t* keys, const rns_params_t* params) {
  poly_free(&keys->s_hat);
  rns_free(&keys->s, params);
  for (uint32_t i = 0; i < params->l; i++) {
    rlwe_free(&keys->r[i], params);
  }
  free(keys->r);
  uint32_t log2_n = floor_log2(params->n) - 1;
  for (uint32_t i = 0; i < log2_n; i++) {
    for (uint32_t j = 0; j < params->l; j++) {
      rlwe_free(&keys->h[i][j], params);
    }
    free(keys->h[i]);
  }
  free(keys->h);
}

void rlwe_gen_secret_key(keys_t* keys, const rns_params_t* params) {
  poly_ternary_sampler(keys->s_hat, params->n);
  for (uint32_t i = 0; i < params->l; i++) {
    memcpy(keys->s[i], keys->s_hat, sizeof(int64_t) * params->n);
    poly_ntt(keys->s[i], params->w[i], params->q[i], params->qRq[i],
             params->qinv[i], params->n);
  }
}

void rlwe_gen_relinearization_key(keys_t* keys, const rns_params_t* params) {
  poly_t s2, e, me;
  poly_new(&s2, params->n);
  poly_new(&e, params->n);
  poly_new(&me, params->n);
  // Compute S^2 in the first (biggest) NTT domain
  poly_pointwise_mul(s2, keys->s[0], keys->s[0], params->q[0], params->qinv[0],
                     params->n);
  poly_intt(s2, params->w_1[0], params->q[0], params->qRq[0], params->qinv[0],
            params->n);
  // Removes the INTT R^2/n factor
  poly_mulx(s2, params->qR2n[0], params->q[0], params->qinv[0], params->n);
  rlwe_key_switch_encrypt(keys->r, s2, me, e, keys, params);
  poly_free(&s2);
  poly_free(&e);
  poly_free(&me);
}

// Trade-off memory/performance for inner-product:
// (log2_n - 1) keys + (log2_n - 1) rotations for reduction
// log2_n - 1 because we work with n/2 complex vector slots
// Our rotation keys are 5^(n/4), 5^(n/4+n/8), ..., 5^(n/2-1)
// We use these keys to reduce the memory usage during reductions
void rlwe_gen_rotation_key(keys_t* keys, const rns_params_t* params) {
  poly_t s2, e, me;
  poly_new(&s2, params->n);
  poly_new(&e, params->n);
  poly_new(&me, params->n);
  uint32_t log2_n2 = floor_log2(params->n) - 1;
  uint32_t exp = 0;
  for (uint32_t i = 0; i < log2_n2; i++) {
    exp += params->n >> (i + 2);
    uint32_t g = mod_pow2(5, exp, 2 * params->n);
    poly_permute(s2, keys->s_hat, g, params->n);
    rlwe_key_switch_encrypt(keys->h[i], s2, me, e, keys, params);
  }
  poly_free(&s2);
  poly_free(&e);
  poly_free(&me);
}

// RLWE encryption step (one of the l polynomials)
static void rlwe_encrypt_step(poly_t a, poly_t b, poly_t m2, const poly_t m,
                              const poly_t s, const poly_t w, uint64_t q,
                              int64_t qRq, int64_t qinv, uint32_t n) {
  // B = -(-A)S + DeltaM + E
  poly_mod(m2, m, q, qRq, n);
  poly_ntt(m2, w, q, qRq, qinv, n);
  poly_uniform_sampler(a, q, qRq, n);  // -A
  poly_pointwise_mul(b, a, s, q, qinv, n);
  poly_rsub(b, m2, q, qRq, n);
}

void rlwe_plaintext_encrypt(rlwe_t* ab, const poly_t m, const keys_t* keys,
                            const rns_params_t* params) {
  // DeltaM + E
  poly_t e, me;
  poly_new(&e, params->n);
  poly_new(&me, params->n);
  poly_cbd24_sampler_add(e, m, params->n);
  // B = -(-A)S + DeltaM + E
  for (uint32_t i = 0; i < params->l; i++) {
    rlwe_encrypt_step(ab->a[i], ab->b[i], me, e, keys->s[i], params->w[i],
                      params->q[i], params->qRq[i], params->qinv[i], params->n);
  }
  poly_free(&e);
  poly_free(&me);
}

// RLWE decryption step
static void rlwe_decrypt_step(poly_t m, const poly_t a, const poly_t b,
                              const poly_t s, const poly_t w, uint64_t q,
                              int64_t qRq, int64_t qinv, int64_t qRn,
                              uint32_t n) {
  // B + (-A)S = DeltaM + E
  poly_pointwise_mul(m, a, s, q, qinv, n);
  poly_add(m, b, q, qRq, n);
  poly_intt(m, w, q, qRq, qinv, n);
  poly_mulx_exact(m, qRn, q, qinv, n);
}

void rlwe_plaintext_decrypt(poly_t m, const rlwe_t* ab, const keys_t* keys,
                            const rns_params_t* params) {
  // For decrypting we only need l = 0 (deltaM + E < q_0/2)
  rlwe_decrypt_step(m, ab->a[0], ab->b[0], keys->s[0], params->w_1[0],
                    params->q[0], params->qRq[0], params->qinv[0],
                    params->qRn[0], params->n);
}

// Since this function is normally called several times by other primitives,
// we ask them to allocate e and me
void rlwe_key_switch_encrypt(switch_t key_switch, const poly_t m, poly_t me,
                             poly_t e, const keys_t* keys,
                             const rns_params_t* params) {
  // B = -(-A)S + DeltaM + E
  for (uint32_t k = 0; k < params->l; k++) {
    // E
    poly_cbd24_sampler(e, params->n);
    // Encrypt 0
    for (uint32_t i = 0; i < params->l; i++) {
      if (i == k) continue;
      rlwe_encrypt_step(key_switch[k].a[i], key_switch[k].b[i], me, e,
                        keys->s[i], params->w[i], params->q[i], params->qRq[i],
                        params->qinv[i], params->n);
    }
    // k level (Encrypt s^2)
    poly_add_no_mod(e, m, params->n);
    rlwe_encrypt_step(key_switch[k].a[k], key_switch[k].b[k], me, e, keys->s[k],
                      params->w[k], params->q[k], params->qRq[k],
                      params->qinv[k], params->n);
  }
}

void rlwe_key_switch_decrypt(poly_t m, const switch_t key_switch,
                             const keys_t* keys, const rns_params_t* params) {
  // For decrypting we only need k = 0, l = 0 and any other index
  // Subtract both of them and completely remove the error
  rlwe_decrypt_step(m, key_switch[0].a[0], key_switch[0].b[0], keys->s[0],
                    params->w_1[0], params->q[0], params->qRq[0],
                    params->qinv[0], params->qRn[0], params->n);
  poly_t m2;
  poly_new(&m2, params->n);
  rlwe_decrypt_step(m2, key_switch[0].a[1], key_switch[0].b[1], keys->s[1],
                    params->w_1[1], params->q[1], params->qRq[1],
                    params->qinv[1], params->qRn[1], params->n);
  poly_sub_no_mod(m, m2, params->n);
  poly_free(&m2);
}

void rlwe_key_switch(rlwe_t* ab, const rns_t a, const switch_t r, poly_t tmp,
                     uint32_t level, const rns_params_t* params) {
  for (uint32_t i = 0; i < level; i++) {
    memset(ab->a[i], 0, params->n * sizeof(int64_t));
  }
  for (uint32_t i = 0; i < level; i++) {
    for (uint32_t j = 0; j < level; j++) {
      poly_pointwise_mul(tmp, a[i], r[i].b[j], params->q[j], params->qinv[j],
                         params->n);  // 1
      poly_sub(ab->b[j], tmp, params->q[j], params->qRq[j],
               params->n);  // 1
      poly_pointwise_mul(tmp, a[i], r[i].a[j], params->q[j], params->qinv[j],
                         params->n);  // R
      poly_sub(ab->a[j], tmp, params->q[j], params->qRq[j],
               params->n);  // R
    }
  }
}
