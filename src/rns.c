
#include "rns.h"

#include <stdio.h>
#include <stdlib.h>

#include "mod.h"

void rns_new(rns_t* a, const rns_params_t* params) {
  *a = malloc(sizeof(uint64_t*) * params->l);
  for (uint32_t i = 0; i < params->l; i++) {
    poly_new(&(*a)[i], params->n);
  }
}

void rns_free(rns_t* a, const rns_params_t* params) {
  for (uint32_t i = 0; i < params->l; i++) {
    poly_free(&(*a)[i]);
  }
  free(*a);
}

void rns_mod_switch(rns_t a, poly_t tmp, uint32_t level,
                    const rns_params_t* params) {
  // For ModSwitch we need to come back to the normal domain
  poly_intt(a[level], params->w_1[level], params->q[level], params->qRq[level],
            params->qinv[level], params->n);
  poly_mulx(a[level], params->qRn[level], params->q[level], params->qinv[level],
            params->n);
  for (uint32_t i = 0; i < level; i++) {
    poly_intt(a[i], params->w_1[i], params->q[i], params->qRq[i],
              params->qinv[i], params->n);
    poly_mulx(a[i], params->qRn[i], params->q[i], params->qinv[i], params->n);
    poly_mod(tmp, a[level], params->q[i], params->qRq[i], params->n);
    poly_sub(a[i], tmp, params->q[i], params->qRq[i], params->n);
    poly_mulx(a[i], params->qinvq[level * params->l + i], params->q[i],
              params->qinv[i], params->n);
    poly_ntt(a[i], params->w[i], params->q[i], params->qRq[i], params->qinv[i],
             params->n);
  }
}

static void rns_gen_roots(rns_params_t* params) {
  for (uint32_t i = 0; i < params->l; i++) {
    poly_gen_roots(params->w[i], params->n, (int64_t)params->q[i],
                   params->qR[i], params->qinv[i]);
  }
}

static void rns_gen_inverse_roots(rns_params_t* params) {
  for (uint32_t i = 0; i < params->l; i++) {
    poly_gen_inverse_roots(params->w_1[i], params->w[i][params->n / 2],
                           params->n, (int64_t)params->q[i], params->qinv[i]);
  }
}

static inline uint32_t compute_ring_dimension(uint32_t qtotal_bits) {
  return 16;
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
  params->qR2n = malloc(sizeof(int64_t) * l);
  params->qn = malloc(sizeof(int64_t) * l);
  params->qRn = malloc(sizeof(int64_t) * l);
  params->qR = malloc(sizeof(int64_t) * l);
  params->qR2 = malloc(sizeof(int64_t) * l);
  params->qRq = malloc(sizeof(int64_t) * l);
  params->qinv = malloc(sizeof(int64_t) * l);
  params->qinvq = malloc(sizeof(int64_t) * l * l);
  params->scales = malloc(sizeof(int64_t) * l);
  rns_new(&params->w, params);
  rns_new(&params->w_1, params);
  // 2n: negacyclic
  find_primes(params->q, q0_bits, scale_bits, 2 * n, l);
  // Init data (less roots)
  for (uint32_t i = 0; i < l; i++) {
    // qinv
    params->qinv[i] = (int64_t)mod64_pow_ct(
        params->q[i], (1UL << 63) - 1);  // q-1 = q^(phi-1) mod 2^64
    // qR
    uint64_t R = (uint64_t)(((__uint128_t)1 << 64) - params->q[i]);
    params->qR[i] = (int64_t)(R % params->q[i]);
    if (params->qR[i] > (int64_t)(params->q[i] / 2))
      params->qR[i] -= (int64_t)params->q[i];
    // qR2
    params->qR2[i] =
        (int64_t)(((__int128_t)params->qR[i] * params->qR[i]) % params->q[i]);
    if (params->qR2[i] > (int64_t)(params->q[i] / 2))
      params->qR2[i] -= (int64_t)params->q[i];
    if (params->qR2[i] < (int64_t) - (params->q[i] / 2))
      params->qR2[i] += (int64_t)params->q[i];
    // qn
    params->qn[i] = mod_pow_ct(n, params->q[i] - 2, (int64_t)params->q[i],
                               params->qR[i], params->qinv[i]);  // 1/n
    if (params->qn[i] > (int64_t)(params->q[i] / 2))
      params->qn[i] -= (int64_t)params->q[i];
    if (params->qn[i] < (int64_t) - (params->q[i] / 2))
      params->qn[i] += (int64_t)params->q[i];
    // qRn
    params->qRn[i] = (int64_t)(((__int128_t)params->qR[i] * params->qn[i]) %
                               params->q[i]);  // R/n
    if (params->qRn[i] > (int64_t)(params->q[i] / 2))
      params->qRn[i] -= (int64_t)params->q[i];
    if (params->qRn[i] < (int64_t) - (params->q[i] / 2))
      params->qRn[i] += (int64_t)params->q[i];
    // qR2n
    params->qR2n[i] = (int64_t)(((__int128_t)params->qR[i] * params->qRn[i]) %
                                params->q[i]);  // R^2/n
    if (params->qR2n[i] > (int64_t)(params->q[i] / 2))
      params->qR2n[i] -= (int64_t)params->q[i];
    if (params->qR2n[i] < (int64_t) - (params->q[i] / 2))
      params->qR2n[i] += (int64_t)params->q[i];
    // qRq
    params->qRq[i] =
        (int64_t)((((__uint128_t)1 << 64) + params->q[i] / 2) / params->q[i]);
  }
  // qinvq
  for (uint32_t j = 0; j < l; j++) {
    for (uint32_t i = 0; i < l; i++) {
      if (i == j) {
        params->qinvq[i * l + j] = 0;  // placeholder (never used)
      } else {
        params->qinvq[i * l + j] = mod_pow_ct(
            (int64_t)params->q[i], params->q[j] - 2, (int64_t)params->q[j],
            params->qR[j], params->qinv[j]);  // (q_i)-1
        params->qinvq[i * l + j] =
            (int64_t)(((__int128_t)params->qR[j] * params->qinvq[i * l + j]) %
                      params->q[j]);  // Add R factor
        if (params->qinvq[i * l + j] > (int64_t)(params->q[j] / 2))
          params->qinvq[i * l + j] -= (int64_t)params->q[j];
        if (params->qinvq[i * l + j] < (int64_t) - (params->q[j] / 2))
          params->qinvq[i * l + j] += (int64_t)params->q[j];
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
  free(params->qinv);
  free(params->qinvq);
  free(params->qRq);
  free(params->qR);
  free(params->qR2);
  free(params->q);
  free(params->qR2n);
  free(params->qRn);
  free(params->qn);
  free(params->scales);
  rns_free(&params->w, params);
  rns_free(&params->w_1, params);
}
