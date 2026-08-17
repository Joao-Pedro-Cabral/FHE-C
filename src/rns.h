
#ifndef __RNS_H__
#define __RNS_H__

#include "poly.h"

typedef poly_t* rns_t;  // L polynomials (one per q_i)

typedef struct {
  uint64_t* q;       // l modulus
  int64_t* qinv;     // inverse of the l q's mod 2^64
  int64_t* qR;       // R=2^64 in the l modulis
  int64_t* qR2;      // R^2=2^128 in the l modulis
  int64_t* qRq;      // R/q, round even
  int64_t* qRn;      // R/n, R=2^64 in q modulis
  int64_t* qR2n;     // R^2/n, R=2^64 in q modulis
  uint64_t* scales;  // scales for each moduli
  uint32_t l;        // size of the rns
  uint32_t n;        // size of the polynomials
  rns_t w;           // NTT roots for each ring
  rns_t w_1;         // INTT roots for each ring
} rns_params_t;

/*
  @param a: rns_t to be allocated
  @param params: RNS params
*/
void rns_new(rns_t* a, const rns_params_t* params);

/*
  @param a: rns_t to be set free
  @param params: RNS params
*/
void rns_free(rns_t* a, const rns_params_t* params);

/*
  @param params: RNS params
*/
void rns_params_new(rns_params_t* params, uint32_t q0_bits, uint32_t scale_bits,
                    uint32_t l);

/*
  @param params: RNS params
*/
void rns_params_free(rns_params_t* params);

#endif  // __RNS_H__
