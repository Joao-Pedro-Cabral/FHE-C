
#ifndef __VECTOR_H__
#define __VECTOR_H__

#include <stdint.h>

typedef struct {
  double* rel;  // real part
  double* img;  // imaginary part
} vec;          // N slots

/*
  @param v: complex vector
  @param scale: input Delta
  @param n: size of the vector
*/
void vec_new(vec* v, uint32_t n);

/*
  @param v: complex vector
*/
void vec_free(vec* v);

// TODO: FORWAD hermitian vector
/*
  @param v1: input vector (n/2 elements)
  @param v2: output REVERSED hermitian vector (n elements)
  @param n: number of elements
*/
void vec_hermitian(vec* v2, const vec* v1, uint32_t n);

/*
  @param v2: input REVERSED hermitian vector (n elements)
  @param v1: output vector (n/2 elements)
  @param n: number of elements
*/
void vec_input(vec* v2, const vec* v1, uint32_t n);

/*
  @param w: odd complex powers of unity of 2n
  @param n: size of the underlying ring
*/
void vec_gen_w(vec* w, uint32_t n);

/*
  @param v: input complex vector/output polynomial vector
  @param w: input complex powers of unity
  @param n: size of the polynomial m and the complex vector
*/
void fft(vec* v, const vec* w, uint32_t n);

/*
  @param v: input complex vector/output polynomial vector
  @param w: input complex powers of unity
  @param n: size of the polynomial m and the complex vector
*/
void ifft(vec* v, const vec* w, uint32_t n);

#endif  // __VECTOR_H__
