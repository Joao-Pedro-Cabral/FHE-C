
#ifndef __CKKS_H__
#define __CKKS_H__

#include "rns.h"
#include "vector.h"

/*
  @param v: input n/2 complex vector
  @param scale: input Delta
  @param n: size of the polynomial m (2x v size)
  @param scaled_m: output scaled (delta) polynomial m
*/
void ckks_encode(rns* scaled_m, const vec* v, uint64_t scale, int n);

/*
  @param scaled_m: input scaled (delta) polynomial m
  @param scale: input Delta
  @param n: size of the polynomial m (2x v size)
  @param v: output n/2 complex vector
*/
void ckks_decode(vec* v, const rns* scaled_m, uint64_t scale, int n);

#endif  // __CKKS_H__
