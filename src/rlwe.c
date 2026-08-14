
#include "rlwe.h"

#include <stdlib.h>
#include <string.h>

void rlwe_new(rlwe_t* ab, const rns_params_t* params) {
  rns_new(&ab->b, params);
  rns_new(&ab->a, params);
}

void rlwe_free(rlwe_t* ab, const rns_params_t* params) {
  rns_free(&ab->b, params);
  rns_free(&ab->a, params);
}

void rlwe_keys_new(keys_t* keys, const rns_params_t* params) {
  rns_new(&keys->s, params);
  keys->r = malloc(params->l * sizeof(rlwe_t));
  for (uint32_t i = 0; i < params->l; i++) {
    rlwe_new(&keys->r[i], params);
  }
}

void rlwe_keys_free(keys_t* keys, const rns_params_t* params) {
  rns_free(&keys->s, params);
  for (uint32_t i = 0; i < params->l; i++) {
    rlwe_free(&keys->r[i], params);
  }
  free(keys->r);
}

void rlwe_gen_secret_key(keys_t* keys, const rns_params_t* params) {
  poly_t s_hat;
  poly_new(&s_hat, params->n);
  poly_ternary_sampler(s_hat, params->n);
  for (uint32_t i = 0; i < params->l; i++) {
    memcpy(keys->s[i], s_hat, sizeof(int64_t) * params->n);
    poly_ntt(keys->s[i], params->w[i], params->q[i], params->qRq[i],
             params->qinv[i], params->n);
  }
  poly_free(&s_hat);
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
    // The terms q/q_i and q_i/q cancel in the RNS representation
    // NTT domain -> pointwise operations
    rlwe_encrypt(&keys->r[i], s2, keys, params);
  }
  poly_free(&s2);
}

void rlwe_encrypt(rlwe_t* ab, const poly_t m, const keys_t* keys,
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

void rlwe_decrypt(poly_t m, const rlwe_t* ab, const keys_t* keys,
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
