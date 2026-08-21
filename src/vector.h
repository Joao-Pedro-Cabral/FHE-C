
#ifndef __VECTOR_H__
#define __VECTOR_H__

#include <stdint.h>

typedef struct {
  double* rel;  // real part
  double* img;  // imaginary part
} vec_t;        // m slots

/*
  @param v: complex vector
  @param m: vector slots
*/
void vec_new(vec_t* v, uint32_t m);

/*
  @param v: complex vector
*/
void vec_free(vec_t* v);

/*
  @param w: odd complex powers of unity of 2n
  @param n: size of the underlying ring
*/
void vec_gen_roots(vec_t* w, uint32_t n);

/*
  @param v: input complex vector/output polynomial vector
  @param w: input complex powers of unity
  @param m: vector slots
*/
void fft(vec_t* v, const vec_t* w, uint32_t m);

/*
  @param v: input complex vector/output polynomial vector
  @param w: input complex powers of unity
  @param m: vector slots
*/
void ifft(vec_t* v, const vec_t* w, uint32_t m);

#endif  // __VECTOR_H__
