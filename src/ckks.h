
#ifndef __CKKS_H__
#define __CKKS_H__

#include "rlwe.h"
#include "vector.h"

// Parameters for CKKS usage (should be initialized with ckks_keygen)
typedef struct {
  vec w;           // odd 2n complex roots of unity
  rns_params rns;  // RNS parameters
  keys kkk;        // RLWE context (keys and roots)
  uint64_t scale;  // input scale (Delta)
} ckks_params;

// OpenFHE cyphertext format
typedef struct {
  rlwe ab;                 // RLWE of the cyphertext
  uint32_t n;              // size of the ring
  uint32_t l;              // size of the rns
  uint64_t scale;          // input scale (Delta)
  uint64_t int_scale;      // not used in CKKS (default: 1)
  uint64_t noise_scale;    // not used in CKKS (default: 1)
  uint32_t encoding_type;  // always 4 (CKKS_PACKED_ENCODING)
} cyphertext;

/*
  @param params: CKKS params to be initialized
  @param n: size of the ring
  @param l: size of the rns
  @param q: l primes that decompose ring moduli Q
  @param scale: scale factor
*/
void ckks_keygen(ckks_params* params, uint32_t n, uint32_t l, const uint64_t* q,
                 uint64_t scale);

/*
  @param params: CKKS params to be set free
*/
void ckks_free(ckks_params* params);

/*
  @param v: input n/2 complex vector
  @param params: CKKS params
  @param m: output scaled (delta) polynomial m
*/
void ckks_encode(poly m, const vec* v, const ckks_params* params);

/*
  @param m: input scaled (delta) polynomial m
  @param w: input odd 2n roots of unity
  @param scale: input Delta
  @param n: size of the polynomial m (2x v size)
  @param v: output n/2 complex vector
*/
void ckks_decode(vec* v, const poly m, const ckks_params* params);

/*
  @param m: input polynomial
  @param params: CKKS params
  @param ct: output cyphertext
*/
void ckks_encrypt(cyphertext* ct, const poly m, const ckks_params* params);

/*
  @param ct: input cyphertext
  @param params: CKKS params
  @param m: output polynomial
*/
void ckks_decrypt(poly m, const cyphertext* ct, const ckks_params* params);

/*
  @param params: CKKS params
  @param ct: cyphertext to be allocated
*/
void cyphertext_new(cyphertext* ct, const ckks_params* params);

/*
  @param params: CKKS params
  @param ct: cyphertext to be set free
*/
void cyphertext_free(cyphertext* ct, const ckks_params* params);

#endif  // __CKKS_H__
