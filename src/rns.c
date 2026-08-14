
#include "rns.h"

#include <stdlib.h>

#include "mod.h"

void rns_new(rns* a, const rns_params_t* params) {
  *a = malloc(sizeof(uint64_t*) * params->l);
  for (uint32_t i = 0; i < params->l; i++) {
    poly_new(&(*a)[i], params->n);
  }
}

void rns_free(rns* a, const rns_params_t* params) {
  for (uint32_t i = 0; i < params->l; i++) {
    poly_free(&(*a)[i]);
  }
  free(*a);
}

void rns_gen_roots(rns_params_t* params) {
  for (uint32_t i = 0; i < params->l; i++) {
    poly_gen_roots(params->w[i], params->n, (int64_t)params->q[i],
                   params->qR[i], params->qinv[i]);
  }
}

void rns_gen_inverse_roots(rns_params_t* params) {
  for (uint32_t i = 0; i < params->l; i++) {
    poly_gen_inverse_roots(params->w[i], params->n, (int64_t)params->q[i],
                           params->qinv[i]);
  }
}

void rns_params_new(rns_params_t* params, const uint64_t* q, uint32_t n,
                    uint32_t l) {
  // Allocate Memory
  params->n = n;
  params->l = l;
  params->q = q;
  params->qR2n = malloc(sizeof(int64_t) * l);
  params->qR = malloc(sizeof(int64_t) * l);
  params->qRq = malloc(sizeof(int64_t) * l);
  params->qinv = malloc(sizeof(int64_t) * l);
  rns_new(&params->w, params);
  // Init data (less roots)
  for (uint32_t i = 0; i < l; i++) {
    params->qinv[i] =
        (int64_t)mod64_pow(q[i], (1UL << 63) - 1);  // q-1 = q^(phi-1) mod 2^64
    uint64_t R = (uint64_t)(((__uint128_t)1 << 64) - q[i]);
    params->qR[i] = (int64_t)(R % q[i]);
    if (params->qR[i] > (int64_t)(q[i] / 2)) params->qR[i] -= (int64_t)q[i];
    params->qR2n[i] = (int64_t)mod_pow(n, q[i] - 2, (int64_t)q[i],
                                       params->qR[i], params->qinv[i]);  // 1/n
    params->qR2n[i] =
        (int64_t)((((__int128_t)1 << 64) * params->qR2n[i]) % q[i]);  // R/n
    params->qR2n[i] =
        (int64_t)((((__int128_t)1 << 64) * params->qR2n[i]) % q[i]);  // R^2/n
    if (params->qR2n[i] > (int64_t)(q[i] / 2)) params->qR2n[i] -= (int64_t)q[i];
    if (params->qR2n[i] < (int64_t) - (q[i] / 2))
      params->qR2n[i] += (int64_t)q[i];
    params->qRq[i] = (int64_t)((((__uint128_t)1 << 64) + q[i] / 2) / q[i]);
  }
}

void rns_params_free(rns_params_t* params) {
  free(params->qinv);
  free(params->qRq);
  free(params->qR);
  free(params->qR2n);
  rns_free(&params->w, params);
}
