
#ifndef __VECTOR_H__
#define __VECTOR_H__

#include <stdint.h>

typedef struct {
  int64_t* rel;  // real part
  int64_t* img;  // imaginary part
  int int_part;  // bits to the integer part of the fixed-point representation
} vec;           // N fixed-point slots

/*
  @param v: complex vector
  @param scale: input Delta
  @param n: size of the vector
*/
void vec_new(vec* v, int n, int int_part);

/*
  @param v: complex vector
*/
void vec_free(vec* v);

// TODO: FORWAD hermitian vector
/*
  @param v1: input vector (n/2 elements)
  @param v2: output REVERSED hermitian vector (n elements)
*/
void vec_hermitian(vec* v2, const vec* v1);

/*
  @param v: input complex vector/output polynomial vector
  @param w: input complex powers of unity
  @param n: size of the polynomial m and the complex vector
*/
void fft(vec* v, const vec* w, int n);

/*
  @param v: input complex vector/output polynomial vector
  @param w: input complex powers of unity
  @param n: size of the polynomial m and the complex vector
*/
void ifft(vec* v, const vec* w, int n);

#endif  // __VECTOR_H__
