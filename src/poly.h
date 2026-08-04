
#ifndef __POLY_H__
#define __POLY_H__

#include <stdint.h>

typedef int64_t* poly;  // N integer coefficients

/*
  @param m: polynomial
  @param n: size of the polynomial
*/
void poly_new(poly* m, uint32_t n);

/*
  @param m: polynomial
*/
void poly_free(poly* m);

/*
  @param m2: output polynomial
  @param m: input polynomial
  @param scale: input Delta
  @param int_part: integer part of the original value
  @param n: size of the polynomial
*/
void poly_rescale(poly m2, const poly m, int64_t scale, uint32_t int_part, uint32_t n);

#endif  // __POLY_H__
