
#include "rlwe.h"

#include <stdlib.h>
#include <string.h>

#include "mod.h"
#include "utils.h"
#include "xof.h"

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
  uint32_t lk = params->l + params->k;
  poly_new(&keys->s_hat, n);
  rns_new(&keys->s, params);
  keys->r = malloc(lk * sizeof(rlwe_t));
  for (uint32_t i = 0; i < lk; i++) {
    rlwe_new(&keys->r[i], params);
  }
  rlwe_new(&keys->pk, params);
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
  for (uint32_t i = 0; i < params->l + params->k; i++) {
    rlwe_free(&keys->r[i]);
  }
  free(keys->r);
  for (uint32_t i = 0; i < params->l; i++) {
    rlwe_free(&keys->cg[i]);
  }
  free(keys->cg);
  rlwe_free(&keys->pk);
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

void rlwe_gen_secret_key(keys_t* keys, xof_state_t* xof_state,
                         const rns_params_t* params) {
  uint32_t n = params->n;
  poly_ternary_sampler(keys->s_hat, xof_state, n);
  for (uint32_t i = 0; i < params->l + params->k; i++) {
    memcpy(&keys->s[n * i], keys->s_hat, sizeof(int64_t) * n);
    poly_ntt(&keys->s[n * i], &params->w[n * i], params->q[i], params->R_q[i],
             params->qinvmR[i], n);
  }
}

void rlwe_gen_relinearization_key(keys_t* keys, xof_state_t* xof_state,
                                  const rns_params_t* params) {
  poly_t s2, e, tmp, tmp2;
  uint32_t n = params->n;
  poly_new(&s2, n);
  poly_new(&e, n);
  poly_new(&tmp, n);
  poly_new(&tmp2, n);
  // Compute S^2 in the first (biggest) NTT domain
  poly_pointwise_mul(s2, &keys->s[0], &keys->s[0], params->q[0],
                     params->qinvmR[0], n);
  poly_intt(s2, &params->w_1[0], params->q[0], params->R_q[0],
            params->qinvmR[0], n);
  // Removes the INTT R^2/n factor
  poly_mulx(s2, params->R2_nmq[0], params->q[0], params->qinvmR[0], n);
  // memset(s2, 1, n * sizeof(int64_t));
  rlwe_key_switch_encrypt(keys->r, xof_state, s2, tmp, tmp2, e, keys, params,
                          true, params->l + params->k);
  poly_free(&s2);
  poly_free(&e);
  poly_free(&tmp);
  poly_free(&tmp2);
}

// RLWE encryption step (one of the l polynomials)
static void rlwe_encrypt_step(poly_t a, poly_t b, poly_t m2,
                              xof_state_t* xof_state, const poly_t m,
                              const poly_t s, const poly_t w, uint64_t q,
                              int64_t R_q, int64_t qinvmR, uint32_t n) {
  // B = -AS + DeltaM + E
  poly_mod(m2, m, q, R_q, n);
  poly_ntt(m2, w, q, R_q, qinvmR, n);
  poly_uniform_sampler(a, xof_state, q, R_q, n);  // A
  poly_pointwise_mrsub_exact(b, a, s, m2, q, qinvmR, R_q, n);
}

void rlwe_plaintext_encrypt(rlwe_t* ab, xof_state_t* xof_state, const poly_t m,
                            const keys_t* keys, const rns_params_t* params) {
  // DeltaM + E
  poly_t e, me;
  uint32_t n = params->n;
  poly_new(&e, n);
  poly_new(&me, n);
  uint8_t* buf = malloc(6 * n * sizeof(uint8_t));
  #pragma omp parallel
  {
    poly_cbd24_sampler_add(e, xof_state, buf, m, n);
    // B = -(-A)S + DeltaM + E
    #pragma omp for
    for (uint32_t i = 0; i < params->l; i++) {
      rlwe_encrypt_step(&ab->a[n * i], &ab->b[n * i], me, xof_state, e,
                        &keys->s[n * i], &params->w[n * i], params->q[i],
                        params->R_q[i], params->qinvmR[i], n);
    }
  }
  poly_free(&e);
  poly_free(&me);
  free(buf);
}

// RLWE decryption step
static void rlwe_decrypt_step(poly_t m, const poly_t a, const poly_t b,
                              const poly_t s, const poly_t w, uint64_t q,
                              int64_t R_q, int64_t qinvmR, int64_t R_nmq,
                              uint32_t n) {
  // B + AS = DeltaM + E
  #pragma omp parallel
  {
  poly_pointwise_macc_exact(m, a, s, b, q, qinvmR, R_q, n);
  poly_intt(m, w, q, R_q, qinvmR, n);
  poly_mulx_exact(m, R_nmq, q, qinvmR, n);
  }
}

void rlwe_plaintext_decrypt(poly_t m, const rlwe_t* ab, const keys_t* keys,
                            const rns_params_t* params) {
  // For decrypting we only need l = 0 (deltaM + E < q_0/2)
  rlwe_decrypt_step(m, &ab->a[0], &ab->b[0], &keys->s[0], &params->w_1[0],
                    params->q[0], params->R_q[0], params->qinvmR[0],
                    params->R_nmq[0], params->n);
}

// Since this function is normally called several times by other primitives,
// we ask them to allocate e and me
void rlwe_key_switch_encrypt(switch_t key_switch, xof_state_t* xof_state,
                             const poly_t m, poly_t tmp, poly_t tmp2, poly_t e,
                             const keys_t* keys, const rns_params_t* params,
                             bool mulP, uint32_t l) {
  uint32_t n = params->n;
  // B = -(-A)S + DeltaM + E
  for (uint32_t k = 0; k < l; k++) {
    // E
    poly_cbd24_sampler(e, xof_state, n);
    // Encrypt 0
    for (uint32_t i = 0; i < l; i++) {
      if (i == k) continue;
      rlwe_encrypt_step(&key_switch[k].a[n * i], &key_switch[k].b[n * i], tmp,
                        xof_state, e, &keys->s[n * i], &params->w[n * i],
                        params->q[i], params->R_q[i], params->qinvmR[i], n);
    }
    // k level (Encrypt s^2)
    if (mulP) {
      // Encrypt Ps^2
      if (k < params->l) {
        memcpy(tmp2, m, n * sizeof(int64_t));
        for (uint32_t i = params->l; i < params->l + params->k; i++) {
          poly_mulx(tmp2, (int64_t)params->q[i], params->q[k],
                    params->qinvmR[k], n);
          poly_mulx(tmp2, params->R2mq[k], params->q[k], params->qinvmR[k], n);
        }
        // Encrypt 0
      } else {
        memset(tmp2, 0, n * sizeof(int64_t));
      }
      poly_add_no_mod(e, tmp2, n);
    } else {
      poly_add_no_mod(e, m, n);
    }
    rlwe_encrypt_step(&key_switch[k].a[n * k], &key_switch[k].b[n * k], tmp,
                      xof_state, e, &keys->s[n * k], &params->w[n * k],
                      params->q[k], params->R_q[k], params->qinvmR[k], n);
  }
}

void rlwe_key_switch_decrypt(poly_t m, const switch_t key_switch,
                             const keys_t* keys, const rns_params_t* params) {
  uint32_t n = params->n;
  // For decrypting we only need k = 0, l = 0 and any other index
  // Subtract both of them and completely remove the error
  rlwe_decrypt_step(m, &key_switch[0].a[0], &key_switch[0].b[0], &keys->s[0],
                    &params->w_1[0], params->q[0], params->R_q[0],
                    params->qinvmR[0], params->R_nmq[0], n);
  poly_t m2;
  poly_new(&m2, n);
  rlwe_decrypt_step(m2, &key_switch[0].a[n], &key_switch[0].b[n], &keys->s[n],
                    &params->w_1[n], params->q[1], params->R_q[1],
                    params->qinvmR[1], params->R_nmq[1], n);
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
                         params->qinvmR[j],
                         n);  // 1
      poly_add(&ab->b[n * j], tmp, params->q[j], params->R_q[j],
               n);  // 1
      poly_pointwise_mul(tmp, &a[n * i], &r[i].a[n * j], params->q[j],
                         params->qinvmR[j],
                         n);  // R
      poly_add(&ab->a[n * j], tmp, params->q[j], params->R_q[j],
               n);  // R
    }
  }
}

void rlwe_plaintext_encrypt_server(rlwe_t* ab, xof_state_t* xof_state,
                                   const poly_t m, const rlwe_t pk,
                                   const rns_params_t* params) {
  poly_t u, e1, e2;
  poly_t u_ntt, e1_ntt, e2_ntt;
  uint32_t n = params->n;
  poly_new(&u, n);
  poly_new(&e1, n);
  poly_new(&e2, n);
  poly_new(&u_ntt, n);
  poly_new(&e1_ntt, n);
  poly_new(&e2_ntt, n);
  uint8_t* buf = malloc(6 * n * sizeof(uint8_t));
  poly_ternary_sampler(u, xof_state, n);
  poly_cbd24_sampler_add(e1, xof_state, buf, m, n);
  poly_cbd24_sampler(e2, xof_state, n);
  for (uint32_t i = 0; i < params->l; i++) {
    memcpy(u_ntt, u, sizeof(int64_t) * n);
    poly_ntt(u_ntt, &params->w[n * i], params->q[i], params->R_q[i],
             params->qinvmR[i], n);
    poly_pointwise_mul(&ab->a[n * i], &pk.a[n * i], u_ntt, params->q[i],
                       params->qinvmR[i], n);
    poly_pointwise_mul(&ab->b[n * i], &pk.b[n * i], u_ntt, params->q[i],
                       params->qinvmR[i], n);
    memcpy(e1_ntt, e1, sizeof(int64_t) * n);
    poly_ntt(e1_ntt, &params->w[n * i], params->q[i], params->R_q[i],
             params->qinvmR[i], n);
    poly_add(&ab->b[n * i], e1_ntt, params->q[i], params->R_q[i], n);
    memcpy(e2_ntt, e2, sizeof(int64_t) * n);
    poly_ntt(e2_ntt, &params->w[n * i], params->q[i], params->R_q[i],
             params->qinvmR[i], n);
    poly_add(&ab->a[n * i], e2_ntt, params->q[i], params->R_q[i], n);
  }
  poly_free(&u);
  poly_free(&e1);
  poly_free(&e2);
  poly_free(&u_ntt);
  poly_free(&e1_ntt);
  poly_free(&e2_ntt);
  free(buf);
}
