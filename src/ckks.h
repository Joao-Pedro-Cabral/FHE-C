
#ifndef __CKKS_H__
#define __CKKS_H__

#include "ciphertext.h"
#include "rlwe.h"
#include "vector.h"

// Parameters for CKKS usage (should be initialized with ckks_keygen)
typedef struct {
  vec_t w;           // odd 2n complex roots of unity
  rns_params_t rns;  // RNS parameters
  keys_t keys;       // RLWE context (keys_t and roots)
} ckks_params_t;

/*
  @param params: CKKS params to be initialized
  @param q0_bits: number of bits of q0
  @param scale_bits: number of bits of the scale and the other moduli
  @param l: size of the rns_t
*/
void ckks_keygen(ckks_params_t* params, uint32_t q0_bits, uint32_t scale_bits,
                 uint32_t l);

/*
  @param params: CKKS params
*/
void ckks_relinearization_keygen(ckks_params_t* params);

/*
  @param params: CKKS params
*/
void ckks_conjugation_keygen(ckks_params_t* params);

/*
  @param params: CKKS params
*/
void ckks_public_keygen(ckks_params_t* params);

/*
  @param params: CKKS params
*/
void ckks_inner_product_keygen(ckks_params_t* params);

/*
  @param params: CKKS params
*/
void ckks_linear_transform_keygen(ckks_params_t* params);

/*
  @param params: CKKS params to be set free
*/
void ckks_free(ckks_params_t* params);

/*
  @param v: input n/2 complex vector
  @param params: CKKS params
  @param m: output scaled (delta) polynomial m
*/
void ckks_encode(poly_t m, const vec_t* v, const ckks_params_t* params);

/*
  @param m: input scaled polynomial m
  @param level: current multiplicative depth
  @param params: CKKS params
  @param v: output n/2 complex vector
*/
void ckks_decode(vec_t* v, const poly_t m, uint32_t level,
                 const ckks_params_t* params);

/*
  @param m: input polynomial
  @param params: CKKS params
  @param ct: output ciphertext_t
*/
void ckks_encrypt(ciphertext_t* ct, const poly_t m,
                  const ckks_params_t* params);

/*
  @param ct: input ciphertext_t
  @param params: CKKS params
  @param m: output polynomial
*/
void ckks_decrypt(poly_t m, const ciphertext_t* ct,
                  const ckks_params_t* params);

/*
  @param ct3: output ciphertext_t
  @param ct1: input ciphertext_t
  @param ct2: input ciphertext_t
  @param params: CKKS params
*/
void ckks_ciphertext_ciphertext_addition(ciphertext_t* ct3,
                                         const ciphertext_t* ct1,
                                         const ciphertext_t* ct2,
                                         const ckks_params_t* params);

/*
  @param ct3: output ciphertext_t
  @param ct1: input ciphertext_t
  @param ct2: input ciphertext_t
  @param params: CKKS params
*/
void ckks_ciphertext_ciphertext_multiplication(ciphertext_t* ct3,
                                               const ciphertext_t* ct1,
                                               const ciphertext_t* ct2,
                                               const ckks_params_t* params);

/*
  @param ct2: output ciphertext_t
  @param ct1: input ciphertext_t
  @param m: input plaintext
  @param params: CKKS params
*/
void ckks_ciphertext_plaintext_addition(ciphertext_t* ct2,
                                        const ciphertext_t* ct1, const poly_t m,
                                        const ckks_params_t* params);

/*
  @param ct2: output ciphertext_t
  @param ct1: input ciphertext_t
  @param m: input plaintext
  @param params: CKKS params
*/
void ckks_ciphertext_plaintext_multiplication(ciphertext_t* ct2,
                                              const ciphertext_t* ct1,
                                              const poly_t m,
                                              const ckks_params_t* params);

/*
  @param ct2: output ciphertext_t
  @param ct1: input ciphertext_t
  @param params: CKKS params
*/
void ckks_ciphertext_conjugation(ciphertext_t* ct2, const ciphertext_t* ct1,
                                 const ckks_params_t* params);

/*
  @param ct: output ciphertext_t
  @param m: input plaintext polynomial
  @param params: CKKS params
*/
void ckks_server_encrypt(ciphertext_t* ct, const poly_t m,
                         const ckks_params_t* params);

/*
  @param ct3: output ciphertext_t
  @param ct1: input ciphertext_t
  @param ct2: input ciphertext_t
  @param params: CKKS params
*/
void ckks_ciphertext_inner_product(ciphertext_t* ct3, const ciphertext_t* ct1,
                                   const ciphertext_t* ct2,
                                   const ckks_params_t* params);

/*
  @param ct2: output ciphertext_t
  @param ct1: input ciphertext_t
  @param mat: input linear transform (rotated diagonals)
  @param params: CKKS params
*/
void ckks_ciphertext_linear_transform(ciphertext_t* ct2,
                                      const ciphertext_t* ct1,
                                      const poly_t* mat,
                                      const ckks_params_t* params);

/*
  @param mat2: output matrix plaintext
  @param mat1: input matrix slots
  @param params: CKKS params
*/
void ckks_linear_transform_encode(poly_t* mat2, const vec_t* mat1,
                                  const ckks_params_t* params);

#endif  // __CKKS_H__
