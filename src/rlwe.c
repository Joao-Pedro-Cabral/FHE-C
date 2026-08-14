
#include "rlwe.h"

#include <string.h>

void rlwe_new(rlwe* ab, const rns_params* params) {
  rns_new(&(ab->b), params);
  rns_new(&(ab->a), params);
}

void rlwe_free(rlwe* ab, const rns_params* params) {
  rns_free(&(ab->b), params);
  rns_free(&(ab->a), params);
}

void rlwe_keys_new(keys* kkk, const rns_params* params) {
  rns_new(&(kkk->s), params);
}

void rlwe_keys_free(keys* kkk, const rns_params* params) {
  rns_free(&(kkk->s), params);
}

void rlwe_gen_secret_key(rns s, const rns_params* params) {
  poly s_hat;
  poly_new(&s_hat, params->n);
  poly_ternary_sampler(s_hat, params->n);
  for (uint32_t i = 0; i < params->l; i++) {
    memcpy(s[i], s_hat, sizeof(int64_t) * params->n);
    poly_ntt(s[i], params->w[i], params->q[i], params->qRq[i], params->qinv[i],
             params->n);
  }
  poly_free(&s_hat);
}

void rlwe_encrypt(rlwe* ab, const poly m, const keys* kkk,
                  const rns_params* params) {
  // DeltaM + E
  poly e;
  poly_new(&e, params->n);
  poly_cbd24_sampler_add(e, m, params->n);
  // B = -(-A)S + DeltaM + E
  for (uint32_t i = 0; i < params->l; i++) {
    poly_uniform_sampler(ab->a[i], params->q[i], params->qRq[i],
                         params->n);  // -A
    poly_pointwise_mul(ab->b[i], ab->a[i], kkk->s[i], params->q[i],
                       params->qinv[i], params->n);
    poly_intt(ab->b[i], params->w[i], params->q[i], params->qRq[i],
              params->qinv[i], params->n);
    poly_rsub(ab->b[i], e, params->q[i], params->qRq[i], params->qinv[i],
              params->qR2n[i], params->n);
  }
  poly_free(&e);
}

void rlwe_decrypt(poly m, const rlwe* ab, const keys* kkk,
                  const rns_params* params) {
  // For decrypting we only need l = 0 (deltaM + E < q_0/2)
  // B + (-A)S = DeltaM + E
  poly_pointwise_mul(m, ab->a[0], kkk->s[0], params->q[0], params->qinv[0],
                     params->n);
  poly_intt(m, params->w[0], params->q[0], params->qRq[0], params->qinv[0],
            params->n);
  poly_add(m, ab->b[0], params->q[0], params->qRq[0], params->qinv[0],
           params->qR2n[0], params->n);
}
