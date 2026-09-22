
#include "rlwe.h"

#include <stdlib.h>
#include <string.h>

#include "mod.h"
#include "utils.h"

void rlwe_new(rlwe_t* ab, const rns_params_t* params) {
  rns_new(&ab->b, params);
  rns_new(&ab->a, params);
}

void rlwe_free(rlwe_t* ab) {
  rns_free(&ab->b);
  rns_free(&ab->a);
}

void rlwe_keys_new(keys_t* keys, const rns_params_t* params) {
  uint32_t n = params->n;
  poly_new(&keys->s_hat, n);
  rns_new(&keys->s, params);
  keys->r = malloc(params->l * sizeof(rlwe_t));
  for (uint32_t i = 0; i < params->l; i++) {
    rlwe_new(&keys->r[i], params);
  }
  keys->cg = malloc(params->l * sizeof(rlwe_t));
  for (uint32_t i = 0; i < params->l; i++) {
    rlwe_new(&keys->cg[i], params);
  }
  uint32_t log2_n = floor_log2(n) - 1;
  keys->ip = malloc(log2_n * sizeof(switch_t));
  for (uint32_t i = 0; i < log2_n; i++) {
    keys->ip[i] = malloc(params->l * sizeof(rlwe_t));
    for (uint32_t j = 0; j < params->l; j++) {
      rlwe_new(&keys->ip[i][j], params);
    }
  }
  uint32_t log2_m = floor_log2(n) - 1;
  uint32_t log2_n2 = floor_log2(n / (2 * N1N2_RATIO)) / 2;
  uint32_t log2_n1 = log2_m - log2_n2;
  uint32_t n1 = (1 << log2_n1);
  uint32_t n2 = (1 << log2_n2);
  keys->lt = malloc((n1 + n2) * sizeof(switch_t));
  for (uint32_t i = 0; i < n1 + n2; i++) {
    keys->lt[i] = malloc(params->l * sizeof(rlwe_t));
    for (uint32_t j = 0; j < params->l; j++) {
      rlwe_new(&keys->lt[i][j], params);
    }
  }
}

void rlwe_keys_free(keys_t* keys, const rns_params_t* params) {
  uint32_t n = params->n;
  poly_free(&keys->s_hat);
  rns_free(&keys->s);
  for (uint32_t i = 0; i < params->l; i++) {
    rlwe_free(&keys->r[i]);
  }
  free(keys->r);
  for (uint32_t i = 0; i < params->l; i++) {
    rlwe_free(&keys->cg[i]);
  }
  free(keys->cg);
  uint32_t log2_n = floor_log2(n) - 1;
  for (uint32_t i = 0; i < log2_n; i++) {
    for (uint32_t j = 0; j < params->l; j++) {
      rlwe_free(&keys->ip[i][j]);
    }
    free(keys->ip[i]);
  }
  free(keys->ip);
  uint32_t log2_m = floor_log2(n) - 1;
  uint32_t log2_n2 = floor_log2(n / (2 * N1N2_RATIO)) / 2;
  uint32_t log2_n1 = log2_m - log2_n2;
  uint32_t n1 = (1 << log2_n1);
  uint32_t n2 = (1 << log2_n2);
  for (uint32_t i = 0; i < n1 + n2; i++) {
    for (uint32_t j = 0; j < params->l; j++) {
      rlwe_free(&keys->lt[i][j]);
    }
    free(keys->lt[i]);
  }
  free(keys->lt);
}

void rlwe_gen_secret_key(keys_t* keys, const rns_params_t* params) {
  uint32_t n = params->n;
  poly_ternary_sampler(keys->s_hat, n);
  for (uint32_t i = 0; i < params->l; i++) {
    memcpy(&keys->s[n * i], keys->s_hat, sizeof(int64_t) * n);
    poly_ntt(&keys->s[n * i], &params->w[n * i], params->q[i], params->qRq[i],
             params->qinv[i], n);
  }
}

void rlwe_gen_relinearization_key(keys_t* keys, const rns_params_t* params) {
  poly_t s2, e, me;
  uint32_t n = params->n;
  poly_new(&s2, n);
  poly_new(&e, n);
  poly_new(&me, n);
  // Compute S^2 in the first (biggest) NTT domain
  poly_pointwise_mul(s2, &keys->s[0], &keys->s[0], params->q[0],
                     params->qinv[0], n);
  poly_intt(s2, &params->w_1[0], params->q[0], params->qRq[0], params->qinv[0],
            n);
  // Removes the INTT R^2/n factor
  poly_mulx(s2, params->qR2n[0], params->q[0], params->qinv[0], n);
  rlwe_key_switch_encrypt(keys->r, s2, me, e, keys, params);
  poly_free(&s2);
  poly_free(&e);
  poly_free(&me);
}

// RLWE encryption step (one of the l polynomials)
static void rlwe_encrypt_step(poly_t a, poly_t b, poly_t m2, const poly_t m,
                              const poly_t s, const poly_t w, uint64_t q,
                              int64_t qRq, int64_t qinv, uint32_t n) {
  // B = -AS + DeltaM + E
  poly_mod(m2, m, q, qRq, n);
  poly_ntt(m2, w, q, qRq, qinv, n);
  poly_uniform_sampler(a, q, qRq, n);  // A
  poly_pointwise_mul(b, a, s, q, qinv, n);
  poly_rsub(b, m2, q, qRq, n);
}

void rlwe_plaintext_encrypt(rlwe_t* ab, const poly_t m, const keys_t* keys,
                            const rns_params_t* params) {
  // DeltaM + E
  poly_t e, me;
  uint32_t n = params->n;
  poly_new(&e, n);
  poly_new(&me, n);
  poly_cbd24_sampler_add(e, m, n);
  // B = -(-A)S + DeltaM + E
  for (uint32_t i = 0; i < params->l; i++) {
    rlwe_encrypt_step(&ab->a[n * i], &ab->b[n * i], me, e, &keys->s[n * i],
                      &params->w[n * i], params->q[i], params->qRq[i],
                      params->qinv[i], n);
  }
  poly_free(&e);
  poly_free(&me);
}

// RLWE decryption step
static void rlwe_decrypt_step(poly_t m, const poly_t a, const poly_t b,
                              const poly_t s, const poly_t w, uint64_t q,
                              int64_t qRq, int64_t qinv, int64_t qRn,
                              uint32_t n) {
  // B + AS = DeltaM + E
  poly_pointwise_mul(m, a, s, q, qinv, n);
  poly_add(m, b, q, qRq, n);
  poly_intt(m, w, q, qRq, qinv, n);
  poly_mulx_exact(m, qRn, q, qinv, n);
}

void rlwe_plaintext_decrypt(poly_t m, const rlwe_t* ab, const keys_t* keys,
                            const rns_params_t* params) {
  // For decrypting we only need l = 0 (deltaM + E < q_0/2)
  rlwe_decrypt_step(m, &ab->a[0], &ab->b[0], &keys->s[0], &params->w_1[0],
                    params->q[0], params->qRq[0], params->qinv[0],
                    params->qRn[0], params->n);
}

// Since this function is normally called several times by other primitives,
// we ask them to allocate e and me
void rlwe_key_switch_encrypt(switch_t key_switch, const poly_t m, poly_t me,
                             poly_t e, const keys_t* keys,
                             const rns_params_t* params) {
  uint32_t n = params->n;
  // B = -(-A)S + DeltaM + E
  for (uint32_t k = 0; k < params->l; k++) {
    // E
    poly_cbd24_sampler(e, n);
    // Encrypt 0
    for (uint32_t i = 0; i < params->l; i++) {
      if (i == k) continue;
      rlwe_encrypt_step(&key_switch[k].a[n * i], &key_switch[k].b[n * i], me, e,
                        &keys->s[n * i], &params->w[n * i], params->q[i],
                        params->qRq[i], params->qinv[i], n);
    }
    // k level (Encrypt s^2)
    poly_add_no_mod(e, m, n);
    rlwe_encrypt_step(&key_switch[k].a[n * k], &key_switch[k].b[n * k], me, e,
                      &keys->s[n * k], &params->w[n * k], params->q[k],
                      params->qRq[k], params->qinv[k], n);
  }
}

void rlwe_key_switch_decrypt(poly_t m, const switch_t key_switch,
                             const keys_t* keys, const rns_params_t* params) {
  uint32_t n = params->n;
  // For decrypting we only need k = 0, l = 0 and any other index
  // Subtract both of them and completely remove the error
  rlwe_decrypt_step(m, &key_switch[0].a[0], &key_switch[0].b[0], &keys->s[0],
                    &params->w_1[0], params->q[0], params->qRq[0],
                    params->qinv[0], params->qRn[0], n);
  poly_t m2;
  poly_new(&m2, n);
  rlwe_decrypt_step(m2, &key_switch[0].a[n], &key_switch[0].b[n], &keys->s[n],
                    &params->w_1[n], params->q[1], params->qRq[1],
                    params->qinv[1], params->qRn[1], n);
  poly_sub_no_mod(m, m2, n);
  poly_free(&m2);
}

void rlwe_key_switch(rlwe_t* ab, const rns_t a, const switch_t r, poly_t tmp,
                     uint32_t level, const rns_params_t* params) {
  uint32_t n = params->n;
  for (uint32_t i = 0; i < level; i++) {
    memset(&ab->a[n * i], 0, n * sizeof(int64_t));
  }
  for (uint32_t i = 0; i < level; i++) {
    for (uint32_t j = 0; j < level; j++) {
      poly_pointwise_mul(tmp, &a[n * i], &r[i].b[n * j], params->q[j],
                         params->qinv[j],
                         n);  // 1
      poly_add(&ab->b[n * j], tmp, params->q[j], params->qRq[j],
               n);  // 1
      poly_pointwise_mul(tmp, &a[n * i], &r[i].a[n * j], params->q[j],
                         params->qinv[j],
                         n);  // R
      poly_add(&ab->a[n * j], tmp, params->q[j], params->qRq[j],
               n);  // R
    }
  }
}
