
#ifndef __POLY_H__
#define __POLY_H__

#include <stdint.h>

#include "vector.h"

typedef int64_t* poly_t;  // N integer coefficients

#define SEED_BYTES (128 / 8)

/*
  @param m: polynomial
  @param n: size of the polynomial
*/
void poly_new(poly_t* m, uint32_t n);

/*
  @param m: polynomial
*/
void poly_free(poly_t* m);

// Obs: I know that what really performs the vector <-> poly_t operations is the
// FFT, but I name the functions like that because of how the C code works

/*
  @param v: input vector
  @param m: output polynomial
  @param scale: input Delta
  @param n: size of the polynomial
*/
void vector_to_poly(poly_t m, const vec_t* v, uint64_t scale, uint32_t n);

/*
  @param v: output vector
  @param m: input polynomial
  @param scale: input Delta
  @param n: size of the polynomial
*/
void poly_to_vector(vec_t* v, const poly_t m, uint64_t scale, uint32_t n);

/*
  @param s: output secret polynomial
  @param n: size of the polynomial
*/
void poly_ternary_sampler(poly_t s, uint32_t n);

/*
  @param a: output public polynomial
  @param q: input modulus
  @param qRq: R/q, round even
  @param n: size of the polynomials
*/
void poly_uniform_sampler(poly_t a, uint64_t q, int64_t qRq, uint32_t n);

/*
  @param e: output error polynomial
  @param m: input message polynomial
  @param n: size of the polynomial
  // To mimic a gaussian with std deviation 3.2 with a cbd
  // The minimal cbd parameter is 21
  // However, for efficiency reasons I do a 24 cbd
*/
void poly_cbd24_sampler_add(poly_t e, const poly_t m, uint32_t n);

/*
  @param e: output error polynomial
  @param n: size of the polynomial
  // To mimic a gaussian with std deviation 3.2 with a cbd
  // The minimal cbd parameter is 21
  // However, for efficiency reasons I do a 24 cbd
*/
void poly_cbd24_sampler(poly_t e, uint32_t n);

/*
  @param a: input/output polynomial
  @param b: input polynomial
  @param n: size of the polynomials
*/
void poly_add_no_mod(poly_t a, const poly_t b, uint32_t n);

/*
  @param a: input/output polynomial
  @param b: input polynomial
  @param n: size of the polynomials
*/
void poly_sub_no_mod(poly_t a, const poly_t b, uint32_t n);

/*
  @param a: input/output polynomial
  @param b: input polynomial
  @param q: input moduli
  @param qRq: R/q, round even
  @param n: size of the polynomials
*/
void poly_add(poly_t a, const poly_t b, uint64_t q, int64_t qRq, uint32_t n);

/*
  @param c: output polynomial
  @param a: input polynomial
  @param b: input polynomial
  @param q: input moduli
  @param qRq: R/q, round even
  @param n: size of the polynomials
*/
void poly_add2(poly_t c, const poly_t a, const poly_t b, uint64_t q,
               int64_t qRq, uint32_t n);

/*
  @param a: input/output polynomial
  @param b: input polynomial
  @param q: input moduli
  @param qRq: R/q, round even
  @param n: size of the polynomials
*/
void poly_sub(poly_t a, const poly_t b, uint64_t q, int64_t qRq, uint32_t n);

/*
  @param a: input/output polynomial
  @param b: input polynomial
  @param q: input moduli
  @param qRq: R/q, round even
  @param n: size of the polynomials
*/
void poly_rsub(poly_t a, const poly_t b, uint64_t q, int64_t qRq, uint32_t n);

/*
  @param a: input/output polynomial
  @param b: input polynomial
  @param h: permute from X to X^h
  @param n: size of the polynomials
*/
void poly_permute(poly_t a, const poly_t b, uint32_t h, uint32_t n);

/*
  @param a: input/output polynomial
  @param b: input polynomial
  @param q: input moduli
  @param qinv: input q^(-1) mod 2^64
  @param n: size of the polynomials
*/
void poly_pointwise_mul(poly_t c, const poly_t a, const poly_t b, uint64_t q,
                        int64_t qinv, uint32_t n);

/*
  @param a: input/output polynomial
  @param x: input multiplicative constant
  @param q: input moduli
  @param qinv: input q^(-1) mod 2^64
  @param n: size of the polynomials
*/
void poly_mulx(poly_t a, int64_t x, uint64_t q, int64_t qinv, uint32_t n);

/*
  @param a: input/output polynomial
  @param x: input multiplicative constant
  @param q: input moduli
  @param qinv: input q^(-1) mod 2^64
  @param n: size of the polynomials
*/
void poly_mulx_exact(poly_t a, int64_t x, uint64_t q, int64_t qinv, uint32_t n);

/*
  @param a: input polynomial
  @param b: output polynomial
  @param q: input moduli
  @param qRq: R/q, round even
  @param n: size of the polynomials
*/
void poly_mod(poly_t b, const poly_t a, uint64_t q, int64_t qRq, uint32_t n);

/*
  @param a: input/output polynomial
  @param q: input moduli
  @param qRq: input R=2^64 mod q
  @param qinv: input q^(-1) mod 2^64
  @param n: size of the polynomial
*/
void poly_ntt(poly_t a, const poly_t roots, uint64_t q, int64_t qRq,
              int64_t qinv, uint32_t n);

/*
  @param a: input/output polynomial
  @param q: input moduli
  @param qRq: input R=2^64 mod q
  @param qinv: input q^(-1) mod 2^64
  @param n: size of the polynomial
*/
void poly_intt(poly_t a, const poly_t roots, uint64_t q, int64_t qRq,
               int64_t qinv, uint32_t n);

/*
  @param q: input moduli
  @param qR: input R=2^64 mod q
  @param qinv: input q^(-1) mod 2^64
  @param n: size of the ring
  @param roots: output primitive roots the moduli
*/
void poly_gen_roots(poly_t roots, uint32_t n, int64_t q, int64_t qR,
                    int64_t qinv);

/*
  Suppose that roots was initialized by a previous call to poly_gen_roots
  @param q: input moduli
  @param qinv: input q^(-1) mod 2^64
  @param n: size of the ring
  @param root0: roots[n/2] from poly_gen_roots
  @param roots: output primitive roots the moduli
*/
void poly_gen_inverse_roots(poly_t roots, int64_t root0, uint32_t n, int64_t q,
                            int64_t qinv);

#endif  // __POLY_H__
