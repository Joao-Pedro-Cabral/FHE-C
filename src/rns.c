
#include "rns.h"

#include <stdio.h>
#include <stdlib.h>

#include "mod.h"

void rns_new(rns_t* a, const rns_params_t* params) {
  *a = malloc(sizeof(uint64_t*) * params->l * params->n);
}

void rns_free(rns_t* a) { free(*a); }

void rns_mod_switch(rns_t a, poly_t tmp, uint32_t level,
                    const rns_params_t* params) {
  // For ModSwitch we need to come back to the normal domain
  uint32_t n = params->n;
  // FastBConv mod q_i
  poly_intt(&a[n * level], &params->w_1[n * level], params->q[level],
            params->R_q[level], params->qinvmR[level], n);
  poly_mulx(&a[n * level], params->R_nmq[level], params->q[level],
            params->qinvmR[level], n);
  for (uint32_t i = 0; i < level; i++) {
    poly_intt(&a[n * i], &params->w_1[n * i], params->q[i], params->R_q[i],
              params->qinvmR[i], n);
    poly_mulx(&a[n * i], params->R_nmq[i], params->q[i], params->qinvmR[i], n);
    // FastBConv mod b_j
    poly_mod(tmp, &a[n * level], params->q[i], params->R_q[i], n);
    //
    poly_sub(&a[n * i], tmp, params->q[i], params->R_q[i], n);
    poly_mulx(&a[n * i], params->Rqinvmq[level * params->l + i], params->q[i],
              params->qinvmR[i], n);
    poly_ntt(&a[n * i], &params->w[n * i], params->q[i], params->R_q[i],
             params->qinvmR[i], n);
  }
}

void rns_mod_switch_p(rns_t a, poly_t tmp1, poly_t tmp2, uint32_t level,
                      const rns_params_t* params) {
  // For ModSwitch we need to come back to the normal domain
  uint32_t n = params->n;
  // FastBConv mod q_i
  uint32_t l = params->l;
  poly_intt(&a[n * l], &params->w_1[n * l], params->q[l], params->R_q[l],
            params->qinvmR[l], n);
  poly_mulx(&a[n * l], params->R_nmq[l], params->q[l], params->qinvmR[l], n);
  poly_mulx(&a[n * l], params->Rqinvmq[(l + 1) * (l + 2) + l], params->q[l],
            params->qinvmR[l], n);
  poly_intt(&a[n * (l + 1)], &params->w_1[n * (l + 1)], params->q[(l + 1)],
            params->R_q[(l + 1)], params->qinvmR[(l + 1)], n);
  poly_mulx(&a[n * (l + 1)], params->R_nmq[(l + 1)], params->q[(l + 1)],
            params->qinvmR[(l + 1)], n);
  poly_mulx(&a[n * (l + 1)], params->Rqinvmq[l * (l + 2) + (l + 1)],
            params->q[(l + 1)], params->qinvmR[(l + 1)], n);
  for (uint32_t i = 0; i < level; i++) {
    poly_intt(&a[n * i], &params->w_1[n * i], params->q[i], params->R_q[i],
              params->qinvmR[i], n);
    poly_mulx(&a[n * i], params->R_nmq[i], params->q[i], params->qinvmR[i], n);
    // FastBConv mod b_j
    poly_mulx2(tmp1, &a[n * l], params->q[l + 1], params->q[i], params->qinvmR[i],
               n);
    poly_mulx2(tmp2, &a[n * (l + 1)], params->q[l], params->q[i],
               params->qinvmR[i], n);
    poly_add(tmp1, tmp2, params->q[i], params->R_q[i], n);
    poly_mulx(tmp1, params->R2mq[i], params->q[i], params->qinvmR[i], n);
    //
    poly_sub(&a[n * i], tmp1, params->q[i], params->R_q[i], n);
    poly_mulx(&a[n * i], params->Rqinvmq[level * (l + 2) + i], params->q[i],
              params->qinvmR[i], n);
    poly_ntt(&a[n * i], &params->w[n * i], params->q[i], params->R_q[i],
             params->qinvmR[i], n);
  }
}

static void rns_gen_roots(rns_params_t* params) {
  uint32_t n = params->n;
  for (uint32_t i = 0; i < params->l; i++) {
    poly_gen_roots(&params->w[n * i], n, (int64_t)params->q[i], params->Rmq[i],
                   params->qinvmR[i]);
  }
}

static void rns_gen_inverse_roots(rns_params_t* params) {
  uint32_t n = params->n;
  for (uint32_t i = 0; i < params->l; i++) {
    poly_gen_inverse_roots(&params->w_1[n * i], params->w[n * i + n / 2], n,
                           (int64_t)params->q[i], params->qinvmR[i]);
  }
}

static inline uint32_t compute_ring_dimension(uint32_t qtotal_bits) {
  return 512;
}

void rns_params_new(rns_params_t* params, uint32_t q0_bits, uint32_t scale_bits,
                    uint32_t l) {
  // Allocate Memory
  uint32_t qtotal_bits = q0_bits + l * scale_bits;
  l++;  // q0 prime + l for multiplicative depth
  uint32_t n = compute_ring_dimension(qtotal_bits);
  params->n = n;
  params->l = l;
  params->q = malloc(sizeof(int64_t) * l);
  params->R2_nmq = malloc(sizeof(int64_t) * l);
  params->_nmq = malloc(sizeof(int64_t) * l);
  params->R_nmq = malloc(sizeof(int64_t) * l);
  params->Rmq = malloc(sizeof(int64_t) * l);
  params->R2mq = malloc(sizeof(int64_t) * l);
  params->R_q = malloc(sizeof(int64_t) * l);
  params->qinvmR = malloc(sizeof(int64_t) * l);
  params->Rqinvmq = malloc(sizeof(int64_t) * l * l);
  params->scales = malloc(sizeof(int64_t) * l);
  rns_new(&params->w, params);
  rns_new(&params->w_1, params);
  // 2n: negacyclic
  find_primes(params->q, q0_bits, scale_bits, 2 * n, l);
  // Init data (less roots)
  for (uint32_t i = 0; i < l; i++) {
    // qinvmR
    params->qinvmR[i] = (int64_t)mod64_pow_ct(
        params->q[i], (1UL << 63) - 1);  // q-1 = q^(phi-1) mod 2^64
    // Rmq
    uint64_t R = (uint64_t)(((__uint128_t)1 << 64) - params->q[i]);
    params->Rmq[i] = (int64_t)(R % params->q[i]);
    if (params->Rmq[i] > (int64_t)(params->q[i] / 2))
      params->Rmq[i] -= (int64_t)params->q[i];
    // R2mq
    params->R2mq[i] =
        (int64_t)(((__int128_t)params->Rmq[i] * params->Rmq[i]) % params->q[i]);
    if (params->R2mq[i] > (int64_t)(params->q[i] / 2))
      params->R2mq[i] -= (int64_t)params->q[i];
    if (params->R2mq[i] < (int64_t) - (params->q[i] / 2))
      params->R2mq[i] += (int64_t)params->q[i];
    // _nmq
    params->_nmq[i] = mod_pow_ct(n, params->q[i] - 2, (int64_t)params->q[i],
                               params->Rmq[i], params->qinvmR[i]);  // 1/n
    if (params->_nmq[i] > (int64_t)(params->q[i] / 2))
      params->_nmq[i] -= (int64_t)params->q[i];
    if (params->_nmq[i] < (int64_t) - (params->q[i] / 2))
      params->_nmq[i] += (int64_t)params->q[i];
    // R_nmq
    params->R_nmq[i] = (int64_t)(((__int128_t)params->Rmq[i] * params->_nmq[i]) %
                               params->q[i]);  // R/n
    if (params->R_nmq[i] > (int64_t)(params->q[i] / 2))
      params->R_nmq[i] -= (int64_t)params->q[i];
    if (params->R_nmq[i] < (int64_t) - (params->q[i] / 2))
      params->R_nmq[i] += (int64_t)params->q[i];
    // R2_nmq
    params->R2_nmq[i] = (int64_t)(((__int128_t)params->Rmq[i] * params->R_nmq[i]) %
                                params->q[i]);  // R^2/n
    if (params->R2_nmq[i] > (int64_t)(params->q[i] / 2))
      params->R2_nmq[i] -= (int64_t)params->q[i];
    if (params->R2_nmq[i] < (int64_t) - (params->q[i] / 2))
      params->R2_nmq[i] += (int64_t)params->q[i];
    // R_q
    params->R_q[i] =
        (int64_t)((((__uint128_t)1 << 64) + params->q[i] / 2) / params->q[i]);
  }
  // Rqinvmq
  for (uint32_t j = 0; j < l; j++) {
    for (uint32_t i = 0; i < l; i++) {
      if (i == j) {
        params->Rqinvmq[i * l + j] = 0;  // placeholder (never used)
      } else {
        params->Rqinvmq[i * l + j] = mod_pow_ct(
            (int64_t)params->q[i], params->q[j] - 2, (int64_t)params->q[j],
            params->Rmq[j], params->qinvmR[j]);  // (q_i)-1
        params->Rqinvmq[i * l + j] =
            (int64_t)(((__int128_t)params->Rmq[j] * params->Rqinvmq[i * l + j]) %
                      params->q[j]);  // Add R factor
        if (params->Rqinvmq[i * l + j] > (int64_t)(params->q[j] / 2))
          params->Rqinvmq[i * l + j] -= (int64_t)params->q[j];
        if (params->Rqinvmq[i * l + j] < (int64_t) - (params->q[j] / 2))
          params->Rqinvmq[i * l + j] += (int64_t)params->q[j];
      }
    }
  }
  // scale
  params->scales[l - 1] = (1UL << (scale_bits - 1)) + (1UL << (scale_bits - 2));
  for (uint32_t i = l - 1; i > 0; i--) {
    double scale2 = (double)params->scales[i] * (double)params->scales[i];
    params->scales[i - 1] = (uint64_t)((scale2 + (double)(params->q[i] / 2)) /
                                       (double)params->q[i]);
  }
  // roots
  rns_gen_roots(params);
  rns_gen_inverse_roots(params);
}

void rns_params_free(rns_params_t* params) {
  free(params->qinvmR);
  free(params->Rqinvmq);
  free(params->R_q);
  free(params->Rmq);
  free(params->R2mq);
  free(params->q);
  free(params->R2_nmq);
  free(params->R_nmq);
  free(params->_nmq);
  free(params->scales);
  rns_free(&params->w);
  rns_free(&params->w_1);
}
