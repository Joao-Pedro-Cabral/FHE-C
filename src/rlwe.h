
#ifndef __RLWE_H__
#define __RLWE_H__

#include "rns.h"

typedef struct {
  rns b;  // ciphertext (b = as + deltam + e)
  rns a;  // ciphertext (a = -a)
} rlwe;

// RLWE Keys
typedef struct {
  rns s;  // secret key in l NTT domains
} keys;

/*
  @param params: RNS params
  @param ab: rlwe to be allocated
*/
void rlwe_new(rlwe* ab, const rns_params* params);

/*
  @param params: RNS params
  @param ab: rlwe to be set free
*/
void rlwe_free(rlwe* ab, const rns_params* params);
/*
  @param params: RNS params
  @param kkk: rlwe context to be allocated
*/
void rlwe_keys_new(keys* kkk, const rns_params* params);

/*
  @param params: RNS params
  @param kkk: rlwe context to be set free
*/
void rlwe_keys_free(keys* kkk, const rns_params* params);

/*
  @param s: secret key
  @param params: RNS params
*/
void rlwe_gen_secret_key(rns s, const rns_params* params);

/*
  @param m: input polynomial
  @param kkk: RLWE context
  @param params: RNS params
  @param ab: output encrypted RLWE
*/
void rlwe_encrypt(rlwe* ab, const poly m, const keys* kkk,
                  const rns_params* params);

/*
  @param params: RNS params
  @param kkk: RLWE context
  @param ab: input encrypted RLWE
  @param m: output polynomial
*/
void rlwe_decrypt(poly m, const rlwe* ab, const keys* kkk,
                  const rns_params* params);

#endif  // __RLWE_H__
