
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

static inline uint32_t compute_ring_dimension(uint32_t qtotal_bits) {
  if (qtotal_bits < 27) {
    return 1024;
  } else if (qtotal_bits < 54) {
    return 2048;
  } else if (qtotal_bits < 109) {
    return 4096;
  } else if (qtotal_bits < 218) {
    return 8192;
  } else if (qtotal_bits < 438) {
    return 16384;
  } else if (qtotal_bits < 881) {
    return 32768;
  } else if (qtotal_bits < 1747) {
    return 65536;
  } else if (qtotal_bits < 3523) {
    return 131072;
  }
  printf("Too large q: %d\n!", qtotal_bits);
  return 1024;
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
  params->qR = malloc(sizeof(int64_t) * l);
  params->qRq = malloc(sizeof(int64_t) * l);
  params->qinv = malloc(sizeof(int64_t) * l);
  rns_new(&params->w, params);
  // 2n: negacyclic
  find_primes(params->q, q0_bits, scale_bits, 2 * n, l);
  // Init data (less roots)
  for (uint32_t i = 0; i < l; i++) {
    params->qinv[i] = (int64_t)mod64_pow(
        params->q[i], (1UL << 63) - 1);  // q-1 = q^(phi-1) mod 2^64
    uint64_t R = (uint64_t)(((__uint128_t)1 << 64) - params->q[i]);
    params->qR[i] = (int64_t)(R % params->q[i]);
    if (params->qR[i] > (int64_t)(params->q[i] / 2))
      params->qR[i] -= (int64_t)params->q[i];
    params->qR2n[i] =
        (int64_t)mod_pow(n, params->q[i] - 2, (int64_t)params->q[i],
                         params->qR[i], params->qinv[i]);  // 1/n
    params->qR2n[i] = (int64_t)((((__int128_t)1 << 64) * params->qR2n[i]) %
                                params->q[i]);  // R/n
    params->qR2n[i] = (int64_t)((((__int128_t)1 << 64) * params->qR2n[i]) %
                                params->q[i]);  // R^2/n
    if (params->qR2n[i] > (int64_t)(params->q[i] / 2))
      params->qR2n[i] -= (int64_t)params->q[i];
    if (params->qR2n[i] < (int64_t) - (params->q[i] / 2))
      params->qR2n[i] += (int64_t)params->q[i];
    params->qRq[i] =
        (int64_t)((((__uint128_t)1 << 64) + params->q[i] / 2) / params->q[i]);
  }
}

void rns_params_free(rns_params_t* params) {
  free(params->qinv);
  free(params->qRq);
  free(params->qR);
  free(params->q);
  free(params->qR2n);
  rns_free(&params->w, params);
}
