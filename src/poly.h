
#ifndef __POLY_H__
#define __POLY_H__

#include <stdint.h>

#include "vector.h"

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

// Obs: I know that what really performs the vector <-> poly operations is the
// FFT, but I name the functions like that because of how the C code works

/*
  @param v: input vector
  @param m: output polynomial
  @param scale: input Delta
  @param n: size of the polynomial
*/
void vector_to_poly(poly m, const vec* v, int64_t scale, uint32_t n);

/*
  @param v: output vector
  @param m: input polynomial
  @param scale: input Delta
  @param n: size of the polynomial
*/
void poly_to_vector(vec* v, const poly m, int64_t scale, uint32_t n);

#endif  // __POLY_H__
