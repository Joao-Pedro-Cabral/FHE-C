
#ifndef __RNS_H__
#define __RNS_H__

#include "poly.h"

typedef int64_t* rns_t;  // L polynomials (one per q_i)

// Using dnum = 1, these values are always fine for P
#define P_MOD_BITS 62
#define P_MOD_SIZE 2

// Notation: x_ymz = (x/y) mod z
typedef struct {
  uint64_t* q;       // l+k modulus
  int64_t* qinvmR;   // inverse of the l+k q's mod 2^64
  int64_t* Rqinvmq;  // R(q_i)-1 mod q_j
  int64_t* Rmq;      // R=2^64 in the l+k modulis repeated k times
  int64_t* R2mq;     // R^2=2^128 in the l+k modulis
  int64_t* R_q;      // R/q, round even
  int64_t* _nmq;     // 1/n, R=2^64 in q modulis
  int64_t* R_nmq;    // R/n, R=2^64 in q modulis
  int64_t* R2_nmq;   // R^2/n, R=2^64 in q modulis
  int64_t* zRbq;     // zR's for FastBConv(, b, q): k*l elements
  int64_t* yRbq;     // yR's for FastBConv(, b, q): k*l elements
  int64_t* zRqb;     // zR's for FastBConv(, q, b): l*k*l elements
  int64_t* yRqb;     // yR's for FastBConv(, q, b): l*k*l elements
  uint64_t* scales;  // scales for each moduli
  uint32_t k;        // size of the rns P (ModRaise)
  uint32_t l;        // size of the rns Q (BaseMod)
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
*/
void rns_free(rns_t* a);

/*
  I have z and y here because they could come from 2 different vectors
  @param a: RNS
  @param tmp: auxiliar RNS
  @param l: remaining levels (0-(l-1))
  @param m + k: levels to be removed (m-(m+k-1))
  @param z: fastbconv parameter
  @param y: fastbconv parameter
  @param params: RNS params
*/
void rns_mod_switch(rns_t a, rns_t tmp, uint32_t l, uint32_t m, uint32_t k,
                    const int64_t* z, const int64_t* y,
                    const rns_params_t* params);

/*
  @param a: RNS
  @param level: new level
  @param params: RNS params
*/
void rns_mod_raise(rns_t a, uint32_t level, const rns_params_t* params);

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
