
#ifndef __CKKS_H__
#define __CKKS_H__

#include "rlwe.h"
#include "vector.h"

// Parameters for CKKS usage (should be initialized with ckks_keygen)
typedef struct {
  vec_t w;           // odd 2n complex roots of unity
  rns_params_t rns;  // RNS parameters
  keys_t keys;       // RLWE context (keys_t and roots)
  uint64_t scale;    // input scale (Delta)
} ckks_params_t;

// OpenFHE cyphertext format
typedef struct {
  rlwe_t ab;               // RLWE of the ciphertext_t
  uint32_t n;              // size of the ring
  uint32_t l;              // size of the rns_t
  uint64_t scale;          // input scale (Delta)
  uint64_t int_scale;      // not used in CKKS (default: 1)
  uint64_t noise_scale;    // not used in CKKS (default: 1)
  uint32_t encoding_type;  // always 4 (CKKS_PACKED_ENCODING)
} ciphertext_t;

/*
  @param params: CKKS params to be initialized
  @param q0_bits: number of bits of q0
  @param scale_bits: number of bits of the other moduli (<= q0_bits)
  @param l: size of the rns_t
  @param scale: scale factor
*/
void ckks_keygen(ckks_params_t* params, uint32_t q0_bits, uint32_t scale_bits,
                 uint32_t l, uint64_t scale);

/*
  @param params: CKKS params
*/
void ckks_relinearization_keygen(ckks_params_t* params);

/*
  @param params: CKKS params
*/
void ckks_rotation_keygen(ckks_params_t* params);

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
  @param m: input scaled (delta) polynomial m
  @param w: input odd 2n roots of unity
  @param scale: input Delta
  @param n: size of the polynomial m (2x v size)
  @param v: output n/2 complex vector
*/
void ckks_decode(vec_t* v, const poly_t m, const ckks_params_t* params);

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
  @param params: CKKS params
  @param ct: ciphertext_t to be allocated
*/
void cyphertext_new(ciphertext_t* ct, const ckks_params_t* params);

/*
  @param params: CKKS params
  @param ct: ciphertext_t to be set free
*/
void cyphertext_free(ciphertext_t* ct, const ckks_params_t* params);

#endif  // __CKKS_H__
