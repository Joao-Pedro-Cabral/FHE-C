
#ifndef __AES128CTR_H__
#define __AES128CTR_H__

#ifdef __ARM_NEON__
#include <arm_neon.h>
#else
#include <immintrin.h>
#endif

#include <stddef.h>
#include <stdint.h>

#define AES128CTR_KEYBYTES 16
#define AES128CTR_BLOCKBYTES 64

typedef struct {
#ifdef __ARM_NEON__
  uint8_t rkeys[AES128CTR_KEYBYTES * 11];
  uint8_t n[AES128CTR_KEYBYTES];
#else
  __m128i rkeys[11];
  __m128i n;
#endif
} aes128ctr_ctx;

void aes128ctr_init(aes128ctr_ctx *state, const uint8_t key[AES128CTR_KEYBYTES],
                    uint64_t nonce);

void aes128ctr_squeezeblocks(uint8_t *out, size_t nblocks,
                             aes128ctr_ctx *state);

void aes128ctr_prf(uint8_t *out, size_t outlen,
                   const uint8_t key[AES128CTR_KEYBYTES], uint64_t nonce);

#endif  // __AES128CTR_H__
