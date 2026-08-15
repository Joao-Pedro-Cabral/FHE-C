
#ifndef __RLWE_H__
#define __RLWE_H__

#include "rns.h"

typedef struct {
  rns_t b;  // ciphertext (b = as + deltam + e)
  rns_t a;  // ciphertext (a = -a)
} rlwe_t;

// RLWE Keys
typedef struct {
  rns_t s;    // secret key in l NTT domains
  rlwe_t* r;  // relinearization key (l RLWE's)
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
  @param keys: rlwe_t context to be allocated
*/
void rlwe_keys_new(keys_t* keys, const rns_params_t* params);

/*
  @param params: RNS params
  @param keys: rlwe_t context to be set free
*/
void rlwe_keys_free(keys_t* keys, const rns_params_t* params);

/*
  @param keys: keys
  @param params: RNS params
*/
void rlwe_gen_secret_key(keys_t* keys, const rns_params_t* params);

/*
  @param keys: keys
  @param params: RNS params
*/
void rlwe_gen_relinearization_key(keys_t* keys, const rns_params_t* params);

/*
  @param m: input polynomial
  @param keys: RLWE context
  @param params: RNS params
  @param ab: output encrypted RLWE
*/
void rlwe_plaintext_encrypt(rlwe_t* ab, const poly_t m, const keys_t* keys,
                  const rns_params_t* params);

/*
  @param params: RNS params
  @param keys: RLWE context
  @param ab: input encrypted RLWE
  @param m: output polynomial
*/
void rlwe_plaintext_decrypt(poly_t m, const rlwe_t* ab, const keys_t* keys,
                  const rns_params_t* params);

/*
  @param params: RNS params
  @param keys: RLWE context
  @param ab: input encrypted RLWE
  @param m: output polynomial
*/
void rlwe_relinearization_key_encrypt(rlwe_t* ab, uint32_t index, const poly_t m, const keys_t* keys,
                  const rns_params_t* params);

/*
  @param params: RNS params
  @param keys: RLWE context
  @param ab: input encrypted RLWE
  @param m: output polynomial
*/
void rlwe_relinearization_key_decrypt(poly_t m, uint32_t index, const rlwe_t* ab, const keys_t* keys,
                  const rns_params_t* params);

#endif  // __RLWE_H__
