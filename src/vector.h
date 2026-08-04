
#ifndef __VECTOR_H__
#define __VECTOR_H__

#include <stdint.h>

typedef struct {
  int64_t* rel;  // real part
  int64_t* img;  // imaginary part
  uint32_t int_part;  // bits to the integer part of the fixed-point representation
} vec;           // N fixed-point slots

/*
  @param v: complex vector
  @param scale: input Delta
  @param n: size of the vector
*/
void vec_new(vec* v, uint32_t n, uint32_t int_part);

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
  @param w: input complex powers of unity of n
  @param n: size of the underlying ring
*/
void gen_w(vec* w, uint32_t n);

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
