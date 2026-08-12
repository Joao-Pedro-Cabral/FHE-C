
#ifndef __AES128CTR_H__
#define __AES128CTR_H__

#include <immintrin.h>
#include <stddef.h>
#include <stdint.h>

#define AES128CTR_KEYBYTES 16
#define AES128CTR_BLOCKBYTES 64

typedef struct {
  __m128i rkeys[11];
  __m128i n;
} aes128ctr_ctx;

void aes128ctr_init(aes128ctr_ctx *state, const uint8_t key[AES128CTR_KEYBYTES],
                    uint64_t nonce);

void aes128ctr_squeezeblocks(uint8_t *out, size_t nblocks,
                             aes128ctr_ctx *state);

void aes128ctr_prf(uint8_t *out, size_t outlen,
                   const uint8_t key[AES128CTR_KEYBYTES], uint64_t nonce);

#endif  // __AES128CTR_H__
