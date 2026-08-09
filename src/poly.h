
#ifndef __POLY_H__
#define __POLY_H__

#include <stdint.h>

#include "vector.h"

typedef int64_t* poly;  // N integer coefficients

#define SEED_BYTES (128 / 8)

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
void vector_to_poly(poly m, const vec* v, uint64_t scale, uint32_t n);

/*
  @param v: output vector
  @param m: input polynomial
  @param scale: input Delta
  @param n: size of the polynomial
*/
void poly_to_vector(vec* v, const poly m, uint64_t scale, uint32_t n);

/*
  @param s: output secret polynomial
  @param n: size of the polynomial
*/
void poly_ternary_sampler(poly s, uint32_t n);

/*
  @param a: output public polynomial
  @param q: input modulus
  @param qR: input R=2^64 mod q
  @param qinv: input q^(-1) mod 2^64
  @param n: size of the polynomials
*/
void poly_uniform_sampler(poly a, uint64_t q, int64_t qR, int64_t qinv,
                          uint32_t n);

/*
  @param e: output error polynomial
  @param n: size of the polynomial
  // To mimic a gaussian with std deviation 3.2 with a cbd
  // The minimal cbd parameter is 21
  // However, for efficiency reasons I do a 24 cbd
*/
void poly_cbd24_sampler(poly e, uint32_t n);

/*
  @param a: input/output polynomial
  @param b: input polynomial
  @param q: input moduli
  @param qR: input R=2^64 mod q
  @param qinv: input q^(-1) mod 2^64
  @param n: size of the polynomials
*/
void poly_add(poly a, const poly b, uint64_t q, int64_t qR, int64_t qinv,
              uint32_t n);

/*
  @param a: input/output polynomial
  @param b: input polynomial
  @param n: size of the polynomials
*/
void poly_add_no_mod(poly a, const poly b, uint32_t n);

/*
  @param a: input/output polynomial
  @param q: input moduli
  @param qR: input R=2^64 mod q
  @param qinv: input q^(-1) mod 2^64
  @param n: size of the polynomials
*/
void poly_neg(poly a, uint64_t q, int64_t qR, int64_t qinv, uint32_t n);

/*
  @param a: input/output polynomial
  @param b: input polynomial
  @param q: input moduli
  @param qinv: input q^(-1) mod 2^64
  @param n: size of the polynomials
*/
void poly_pointwise_mul(poly c, const poly a, const poly b, uint64_t q,
                        int64_t qinv, uint32_t n);

/*
  @param a: input/output polynomial
  @param q: input moduli
  @param qR: input R=2^64 mod q
  @param qinv: input q^(-1) mod 2^64
  @param n: size of the polynomial
*/
void poly_ntt(poly a, const poly roots, uint64_t q, int64_t qR, int64_t qinv,
              uint32_t n);

/*
  @param a: input/output polynomial
  @param q: input moduli
  @param qR: input R=2^64 mod q
  @param qinv: input q^(-1) mod 2^64
  @param qn: input R^2/n mod q, R=2^64
  @param n: size of the polynomial
*/
void poly_intt(poly a, const poly roots, uint64_t q, int64_t qR, int64_t qinv,
               int64_t qn, uint32_t n);

/*
  @param q: input moduli
  @param n: size of the ring
  @param roots: output primitive roots the moduli
*/
void poly_gen_roots(poly roots, uint32_t n, int64_t q);

/*
  Suppose that roots was initialized by a previous call to poly_gen_roots
  @param q: input moduli
  @param n: size of the ring
  @param roots: output primitive roots the moduli
*/
void poly_gen_inverse_roots(poly roots, uint32_t n, int64_t q);

#endif  // __POLY_H__
