
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
  poly_t s2;
  poly_new(&s2, params->n);
  // Compute S^2 in the first (biggest) NTT domain
  // The pointwise removes the INTT R^2/n factor
  poly_pointwise_mulx(s2, keys->s[0], keys->s[0], params->qR2n[0], params->q[0],
                      params->qinv[0], params->n);
  poly_intt(s2, params->w[0], params->q[0], params->qRq[0], params->qinv[0],
            params->n);
  for (uint32_t i = 0; i < params->l; i++) {
    rlwe_key_switch_encrypt(&keys->r[i], i, s2, keys, params);
  }
  poly_free(&s2);
}

// Trade-off memory/performance for inner-product:
// (log2_n - 1) keys + (log2_n - 1) rotations for reduction
// log2_n - 1 because we work with complex vector slots
void rlwe_gen_rotation_key(keys_t* keys, const rns_params_t* params) {
  poly_t s2;
  poly_new(&s2, params->n);
  uint32_t log2_n = floor_log2(params->n) - 1;
  uint32_t g = 5;
  uint32_t n2 = 2 * params->n;
  for (uint32_t i = 0; i < log2_n; i++) {
    poly_permute(s2, keys->s_hat, g, params->n);
    for (uint32_t j = 0; j < params->l; j++) {
      rlwe_key_switch_encrypt(&keys->h[i][j], j, s2, keys, params);
    }
    g = (g * g) & (n2 - 1);  // mod 2n
  }
  poly_free(&s2);
}

void rlwe_plaintext_encrypt(rlwe_t* ab, const poly_t m, const keys_t* keys,
                            const rns_params_t* params) {
  // DeltaM + E
  poly_t e;
  poly_new(&e, params->n);
  poly_cbd24_sampler_add(e, m, params->n);
  // B = -(-A)S + DeltaM + E
  for (uint32_t i = 0; i < params->l; i++) {
    poly_uniform_sampler(ab->a[i], params->q[i], params->qRq[i],
                         params->n);  // -A
    poly_pointwise_mul(ab->b[i], ab->a[i], keys->s[i], params->q[i],
                       params->qinv[i], params->n);
    poly_intt(ab->b[i], params->w[i], params->q[i], params->qRq[i],
              params->qinv[i], params->n);
    poly_rsub(ab->b[i], e, params->q[i], params->qRq[i], params->qinv[i],
              params->qR2n[i], params->n);
  }
  poly_free(&e);
}

void rlwe_plaintext_decrypt(poly_t m, const rlwe_t* ab, const keys_t* keys,
                            const rns_params_t* params) {
  // For decrypting we only need l = 0 (deltaM + E < q_0/2)
  // B + (-A)S = DeltaM + E
  poly_pointwise_mul(m, ab->a[0], keys->s[0], params->q[0], params->qinv[0],
                     params->n);
  poly_intt(m, params->w[0], params->q[0], params->qRq[0], params->qinv[0],
            params->n);
  poly_add(m, ab->b[0], params->q[0], params->qRq[0], params->qinv[0],
           params->qR2n[0], params->n);
}

void rlwe_key_switch_encrypt(switch_t key_switch, uint32_t index,
                             const poly_t m, const keys_t* keys,
                             const rns_params_t* params) {
  // DeltaM + E
  poly_t e;
  poly_new(&e, params->n);
  poly_cbd24_sampler(e, params->n);
  // B = -(-A)S + DeltaM + E
  // Encrypt 0
  for (uint32_t i = 0; i < params->l; i++) {
    if (i == index) continue;
    poly_uniform_sampler(key_switch->a[i], params->q[i], params->qRq[i],
                         params->n);  // -A
    poly_pointwise_mul(key_switch->b[i], key_switch->a[i], keys->s[i],
                       params->q[i], params->qinv[i], params->n);
    poly_intt(key_switch->b[i], params->w[i], params->q[i], params->qRq[i],
              params->qinv[i], params->n);
    poly_rsub(key_switch->b[i], e, params->q[i], params->qRq[i],
              params->qinv[i], params->qR2n[i], params->n);
  }
  // index level (Encrypt s^2)
  poly_add_no_mod(e, m, params->n);
  poly_uniform_sampler(key_switch->a[index], params->q[index],
                       params->qRq[index],
                       params->n);  // -A
  poly_pointwise_mul(key_switch->b[index], key_switch->a[index], keys->s[index],
                     params->q[index], params->qinv[index], params->n);
  poly_intt(key_switch->b[index], params->w[index], params->q[index],
            params->qRq[index], params->qinv[index], params->n);
  poly_rsub(key_switch->b[index], e, params->q[index], params->qRq[index],
            params->qinv[index], params->qR2n[index], params->n);
  poly_free(&e);
}

void rlwe_key_switch_decrypt(poly_t m, uint32_t index,
                             const switch_t key_switch, const keys_t* keys,
                             const rns_params_t* params) {
  // For decrypting we only need l = index and any other index
  // Subtract both of them and completely remove the error
  poly_pointwise_mul(m, key_switch->a[index], keys->s[index], params->q[index],
                     params->qinv[index], params->n);
  poly_intt(m, params->w[index], params->q[index], params->qRq[index],
            params->qinv[index], params->n);
  poly_add(m, key_switch->b[index], params->q[index], params->qRq[index],
           params->qinv[index], params->qR2n[index], params->n);
  uint32_t index2 = (index + 1) % params->l;
  poly_t m2;
  poly_new(&m2, params->n);
  poly_pointwise_mul(m2, key_switch->a[index2], keys->s[index2],
                     params->q[index2], params->qinv[index2], params->n);
  poly_intt(m2, params->w[index2], params->q[index2], params->qRq[index2],
            params->qinv[index2], params->n);
  poly_add(m2, key_switch->b[index2], params->q[index2], params->qRq[index2],
           params->qinv[index2], params->qR2n[index2], params->n);
  poly_sub_no_mod(m, m2, params->n);
  poly_free(&m2);
}
