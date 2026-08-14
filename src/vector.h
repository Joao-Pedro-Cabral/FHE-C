
#ifndef __VECTOR_H__
#define __VECTOR_H__

#include <stdint.h>

typedef struct {
  double* rel;  // real part
  double* img;  // imaginary part
} vec_t;        // N slots

/*
  @param v: complex vector
  @param n: size of the vector
*/
void vec_new(vec_t* v, uint32_t n);

/*
  @param v: complex vector
*/
void vec_free(vec_t* v);

// TODO: FORWAD hermitian vector
/*
  @param v1: input vector (n/2 elements)
  @param v2: output REVERSED hermitian vector (n elements)
  @param n: number of elements
*/
void vec_hermitian(vec_t* v2, const vec_t* v1, uint32_t n);

/*
  @param v2: input REVERSED hermitian vector (n elements)
  @param v1: output vector (n/2 elements)
  @param n: number of elements
*/
void vec_input(vec_t* v2, const vec_t* v1, uint32_t n);

/*
  @param w: odd complex powers of unity of 2n
  @param n: size of the underlying ring
*/
void vec_gen_roots(vec_t* w, uint32_t n);

/*
  @param v: input complex vector/output polynomial vector
  @param w: input complex powers of unity
  @param n: size of the polynomial m and the complex vector
*/
void fft(vec_t* v, const vec_t* w, uint32_t n);

/*
  @param v: input complex vector/output polynomial vector
  @param w: input complex powers of unity
  @param n: size of the polynomial m and the complex vector
*/
void ifft(vec_t* v, const vec_t* w, uint32_t n);

#endif  // __VECTOR_H__
