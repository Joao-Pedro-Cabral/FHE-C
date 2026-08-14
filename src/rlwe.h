
#ifndef __RLWE_H__
#define __RLWE_H__

#include "rns.h"

typedef struct {
  rns b;  // ciphertext (b = as + deltam + e)
  rns a;  // ciphertext (a = -a)
} rlwe_t;

// RLWE Keys
typedef struct {
  rns s;  // secret key in l NTT domains
} keys_t;

/*
  @param params: RNS params
  @param ab: rlwe_t to be allocated
*/
void rlwe_new(rlwe_t* ab, const rns_params_t* params);

/*
  @param params: RNS params
  @param ab: rlwe_t to be set free
*/
void rlwe_free(rlwe_t* ab, const rns_params_t* params);
/*
  @param params: RNS params
  @param kkk: rlwe_t context to be allocated
*/
void rlwe_keys_new(keys_t* kkk, const rns_params_t* params);

/*
  @param params: RNS params
  @param kkk: rlwe_t context to be set free
*/
void rlwe_keys_free(keys_t* kkk, const rns_params_t* params);

/*
  @param s: secret key
  @param params: RNS params
*/
void rlwe_gen_secret_key(rns s, const rns_params_t* params);

/*
  @param m: input polynomial
  @param kkk: RLWE context
  @param params: RNS params
  @param ab: output encrypted RLWE
*/
void rlwe_encrypt(rlwe_t* ab, const poly_t m, const keys_t* kkk,
                  const rns_params_t* params);

/*
  @param params: RNS params
  @param kkk: RLWE context
  @param ab: input encrypted RLWE
  @param m: output polynomial
*/
void rlwe_decrypt(poly_t m, const rlwe_t* ab, const keys_t* kkk,
                  const rns_params_t* params);

#endif  // __RLWE_H__
