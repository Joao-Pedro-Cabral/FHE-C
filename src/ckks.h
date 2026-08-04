
#ifndef __CKKS_H__
#define __CKKS_H__

#include "rns.h"
#include "vector.h"

/*
  @param n: size of the polynomial m (2x v size)
  @param w: output odd 2n roots of unity
*/
void ckks_keygen(vec* w, uint32_t n);

/*
  @param v: input n/2 complex vector
  @param w: input odd 2n roots of unity
  @param scale: input Delta
  @param n: size of the polynomial m (2x v size)
  @param scaled_m: output scaled (delta) polynomial m
*/
void ckks_encode(poly scaled_m, const vec* v, const vec* w, int64_t scale, uint32_t n);

/*
  @param scaled_m: input scaled (delta) polynomial m
  @param w: input odd 2n roots of unity
  @param scale: input Delta
  @param n: size of the polynomial m (2x v size)
  @param v: output n/2 complex vector
*/
void ckks_decode(vec* v, const poly scaled_m, const vec* w, int64_t scale, uint32_t n);

#endif  // __CKKS_H__
