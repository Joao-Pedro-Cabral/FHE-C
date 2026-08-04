
#ifndef __POLY_H__
#define __POLY_H__

#include <stdint.h>

typedef int64_t* poly;  // N integer coefficients

/*
  @param m2: output polynomial
  @param m: input polynomial
  @param scale: input Delta
  @param int_part: integer part of the original value
  @param n: size of the vector
*/
void poly_rescale(poly* m2, const poly* m, uint64_t scale, int int_part, int n);

#endif  // __POLY_H__
