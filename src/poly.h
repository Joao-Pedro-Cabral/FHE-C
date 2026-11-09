
#ifndef __POLY_H__
#define __POLY_H__

#include <stdint.h>

#include "vector.h"
#include "xof.h"

typedef int64_t* poly_t;  // N integer coefficients

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
  @param n2: size of the vector slot
  @param n: size of the polynomial
*/
void vector_to_poly(poly_t m, const vec_t* v, uint64_t scale, uint32_t n2,
                    uint32_t n);

/*
  @param v: output vector
  @param m: input polynomial
  @param scale: input Delta
  @param n2: size of the vector slot
  @param n: size of the polynomial
*/
void poly_to_vector(vec_t* v, const poly_t m, uint64_t scale, uint32_t n2,
                    uint32_t n);

/*
  @param s: output secret polynomial
  @param xof_state: XOF state
  @param domain_separator: XOF domain separator
  @param n: size of the polynomial
*/
void poly_ternary_sampler(poly_t s, xof_state_t* xof_state,
                          uint64_t domain_separator, uint8_t* buf, uint32_t n);

/*
  @param a: output public polynomial
  @param xof_state: XOF state
  @param domain_separator: XOF domain separator
  @param q: input modulus
  @param R_q: R/q, round even
  @param n: size of the polynomials
*/
void poly_uniform_sampler(poly_t a, xof_state_t* xof_state,
                          uint64_t domain_separator, uint64_t q, int64_t R_q,
                          uint32_t n);

/*
  @param e: output error polynomial
  @param m: input message polynomial
  @param xof_state: XOF state
  @param domain_separator: XOF domain separator
  @param n: size of the polynomial
  // To mimic a gaussian with std deviation 3.2 with a cbd
  // The minimal cbd parameter is 21
  // However, for efficiency reasons I do a 24 cbd
*/
void poly_cbd24_sampler_add(poly_t e, xof_state_t* xof_state,
                            uint64_t domain_separator, uint8_t* buf,
                            const poly_t m, uint32_t n);

/*
  @param e: output error polynomial
  @param xof_state: XOF state
  @param domain_separator: XOF domain separator
  @param n: size of the polynomial
  // To mimic a gaussian with std deviation 3.2 with a cbd
  // The minimal cbd parameter is 21
  // However, for efficiency reasons I do a 24 cbd
*/
void poly_cbd24_sampler(poly_t e, xof_state_t* xof_state,
                        uint64_t domain_separator, uint32_t n);

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
  @param R_q: R/q, round even
  @param n: size of the polynomials
*/
void poly_add(poly_t a, const poly_t b, uint64_t q, int64_t R_q, uint32_t n);

/*
  @param c: output polynomial
  @param a: input polynomial
  @param b: input polynomial
  @param q: input moduli
  @param R_q: R/q, round even
  @param n: size of the polynomials
*/
void poly_add2(poly_t c, const poly_t a, const poly_t b, uint64_t q,
               int64_t R_q, uint32_t n);

/*
  @param a: input/output polynomial
  @param b: input polynomial
  @param q: input moduli
  @param R_q: R/q, round even
  @param n: size of the polynomials
*/
void poly_sub(poly_t a, const poly_t b, uint64_t q, int64_t R_q, uint32_t n);

/*
  @param a: input/output polynomial
  @param b: input polynomial
  @param q: input moduli
  @param R_q: R/q, round even
  @param n: size of the polynomials
*/
void poly_rsub(poly_t a, const poly_t b, uint64_t q, int64_t R_q, uint32_t n);

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
  @param qinvmR: input q^(-1) mod 2^64
  @param n: size of the polynomials
*/
void poly_pointwise_mul(poly_t c, const poly_t a, const poly_t b, uint64_t q,
                        int64_t qinvmR, uint32_t n);

/*
  @param a: input/output polynomial
  @param b: input polynomial
  @param c: input polynomial
  @param q: input moduli
  @param qinvmR: input q^(-1) mod 2^64
  @param n: size of the polynomials
*/
void poly_pointwise_macc_exact(poly_t d, const poly_t a, const poly_t b,
                               const poly_t c, uint64_t q, int64_t qinvmR,
                               int64_t R_q, uint32_t n);

/*
  @param a: input/output polynomial
  @param b: input polynomial
  @param c: input polynomial
  @param q: input moduli
  @param qinvmR: input q^(-1) mod 2^64
  @param n: size of the polynomials
*/
void poly_pointwise_mrsub_exact(poly_t d, const poly_t a, const poly_t b,
                                const poly_t c, uint64_t q, int64_t qinvmR,
                                int64_t R_q, uint32_t n);

/*
  @param a: input/output polynomial
  @param x: input multiplicative constant
  @param q: input moduli
  @param qinvmR: input q^(-1) mod 2^64
  @param n: size of the polynomials
*/
void poly_mulx(poly_t a, int64_t x, uint64_t q, int64_t qinvmR, uint32_t n);

/*
  @param a: output polynomial
  @param b: input polynomial
  @param x: input multiplicative constant
  @param q: input moduli
  @param qinvmR: input q^(-1) mod 2^64
  @param n: size of the polynomials
*/
void poly_mulx2(poly_t a, const poly_t b, int64_t x, uint64_t q, int64_t qinvmR,
                uint32_t n);

/*
  @param a: input/output polynomial
  @param x: input multiplicative constant
  @param q: input moduli
  @param qinvmR: input q^(-1) mod 2^64
  @param n: size of the polynomials
*/
void poly_mulx_exact(poly_t a, int64_t x, uint64_t q, int64_t qinvmR,
                     uint32_t n);

/*
  @param a: input polynomial
  @param b: output polynomial
  @param q: input moduli
  @param R_q: R/q, round even
  @param n: size of the polynomials
*/
void poly_mod(poly_t b, const poly_t a, uint64_t q, int64_t R_q, uint32_t n);

/*
  @param a: input/output polynomial
  @param b: input polynomial
  @param n: size of the polynomials
*/
void poly_add_no_mod_omp(poly_t a, const poly_t b, uint32_t n);

/*
  @param a: input/output polynomial
  @param b: input polynomial
  @param n: size of the polynomials
*/
void poly_sub_no_mod_omp(poly_t a, const poly_t b, uint32_t n);

/*
  @param a: input/output polynomial
  @param b: input polynomial
  @param q: input moduli
  @param R_q: R/q, round even
  @param n: size of the polynomials
*/
void poly_add_omp(poly_t a, const poly_t b, uint64_t q, int64_t R_q,
                  uint32_t n);

/*
  @param c: output polynomial
  @param a: input polynomial
  @param b: input polynomial
  @param q: input moduli
  @param R_q: R/q, round even
  @param n: size of the polynomials
*/
void poly_add2_omp(poly_t c, const poly_t a, const poly_t b, uint64_t q,
                   int64_t R_q, uint32_t n);

/*
  @param a: input/output polynomial
  @param b: input polynomial
  @param q: input moduli
  @param R_q: R/q, round even
  @param n: size of the polynomials
*/
void poly_sub_omp(poly_t a, const poly_t b, uint64_t q, int64_t R_q,
                  uint32_t n);

/*
  @param a: input/output polynomial
  @param b: input polynomial
  @param q: input moduli
  @param R_q: R/q, round even
  @param n: size of the polynomials
*/
void poly_rsub_omp(poly_t a, const poly_t b, uint64_t q, int64_t R_q,
                   uint32_t n);

/*
  @param a: input/output polynomial
  @param b: input polynomial
  @param h: permute from X to X^h
  @param n: size of the polynomials
*/
void poly_permute_omp(poly_t a, const poly_t b, uint32_t h, uint32_t n);

/*
  @param c: output polynomial
  @param a: input polynomial
  @param b: input polynomial
  @param q: input moduli
  @param qinvmR: input q^(-1) mod 2^64
  @param n: size of the polynomials
*/
void poly_pointwise_mul_omp(poly_t c, const poly_t a, const poly_t b,
                            uint64_t q, int64_t qinvmR, uint32_t n);

/*
  @param d: output polynomial
  @param a: input polynomial
  @param b: input polynomial
  @param c: input polynomial
  @param q: input moduli
  @param qinvmR: input q^(-1) mod 2^64
  @param R_q: R/q, round even
  @param n: size of the polynomials
*/
void poly_pointwise_macc_exact_omp(poly_t d, const poly_t a, const poly_t b,
                                   const poly_t c, uint64_t q, int64_t qinvmR,
                                   int64_t R_q, uint32_t n);

/*
  @param d: output polynomial
  @param a: input polynomial
  @param b: input polynomial
  @param c: input polynomial
  @param q: input moduli
  @param qinvmR: input q^(-1) mod 2^64
  @param R_q: R/q, round even
  @param n: size of the polynomials
*/
void poly_pointwise_mrsub_exact_omp(poly_t d, const poly_t a, const poly_t b,
                                    const poly_t c, uint64_t q, int64_t qinvmR,
                                    int64_t R_q, uint32_t n);

/*
  @param a: input/output polynomial
  @param x: input multiplicative constant
  @param q: input moduli
  @param qinvmR: input q^(-1) mod 2^64
  @param n: size of the polynomials
*/
void poly_mulx_omp(poly_t a, int64_t x, uint64_t q, int64_t qinvmR, uint32_t n);

/*
  @param a: output polynomial
  @param b: input polynomial
  @param x: input multiplicative constant
  @param q: input moduli
  @param qinvmR: input q^(-1) mod 2^64
  @param n: size of the polynomials
*/
void poly_mulx2_omp(poly_t a, const poly_t b, int64_t x, uint64_t q,
                    int64_t qinvmR, uint32_t n);

/*
  @param a: input/output polynomial
  @param x: input multiplicative constant
  @param q: input moduli
  @param qinvmR: input q^(-1) mod 2^64
  @param n: size of the polynomials
*/
void poly_mulx_exact_omp(poly_t a, int64_t x, uint64_t q, int64_t qinvmR,
                         uint32_t n);

/*
  @param b: output polynomial
  @param a: input polynomial
  @param q: input moduli
  @param R_q: R/q, round even
  @param n: size of the polynomials
*/
void poly_mod_omp(poly_t b, const poly_t a, uint64_t q, int64_t R_q,
                  uint32_t n);

/*
  @param a: input/output polynomial
  @param q: input moduli
  @param R_q: input R=2^64 mod q
  @param qinvmR: input q^(-1) mod 2^64
  @param n: size of the polynomial
*/
void poly_ntt(poly_t a, const poly_t roots, uint64_t q, int64_t R_q,
              int64_t qinvmR, uint32_t n);

/*
  @param a: input/output polynomial
  @param q: input moduli
  @param R_q: input R=2^64 mod q
  @param qinvmR: input q^(-1) mod 2^64
  @param n: size of the polynomial
*/
void poly_intt(poly_t a, const poly_t roots, uint64_t q, int64_t R_q,
               int64_t qinvmR, uint32_t n);

/*
  @param q: input moduli
  @param Rmq: input R=2^64 mod q
  @param qinvmR: input q^(-1) mod 2^64
  @param n: size of the ring
  @param roots: output primitive roots the moduli
*/
void poly_gen_roots(poly_t roots, uint32_t n, int64_t q, int64_t Rmq,
                    int64_t qinvmR);

/*
  Suppose that roots was initialized by a previous call to poly_gen_roots
  @param q: input moduli
  @param qinvmR: input q^(-1) mod 2^64
  @param n: size of the ring
  @param root0: roots[n/2] from poly_gen_roots
  @param roots: output primitive roots the moduli
*/
void poly_gen_inverse_roots(poly_t roots, int64_t root0, uint32_t n, int64_t q,
                            int64_t qinvmR);

#endif  // __POLY_H__
