
#ifndef __CKKS_H__
#define __CKKS_H__

#include "rlwe.h"
#include "vector.h"

typedef struct {
  vec* w;         // odd 2n roots of unity
  set q;          // l modulus
  set* eta;       // l sets of roots of unity
  rlwe s;         // secret key in l NTT domains
  rlwe p;         // public key in l NTT domains
  uint32_t l;     // size of the RNS
  uint32_t n;     // size of the polynomials
  int64_t scale;  // input scale (Delta)
} ckks_params;

/*
  @param n: size of the polynomial m (2x v size)
  @param w: output odd 2n roots of unity
*/
void ckks_init(ckks_params* params);

/*
  @param v: input n/2 complex vector
  @param w: input odd 2n roots of unity
  @param scale: input Delta
  @param n: size of the polynomial m (2x v size)
  @param m: output scaled (delta) polynomial m
*/
void ckks_encode(poly m, const vec* v, ckks_params* params);

/*
  @param m: input scaled (delta) polynomial m
  @param w: input odd 2n roots of unity
  @param scale: input Delta
  @param n: size of the polynomial m (2x v size)
  @param v: output n/2 complex vector
*/
void ckks_decode(vec* v, const poly m, ckks_params* params);

#endif  // __CKKS_H__
