
#ifndef __AES128CTR_H__
#define __AES128CTR_H__

#ifdef __NEON_MODE__
#include <arm_neon.h>
#elif __AVX2_MODE__
#include <immintrin.h>
#endif

#include <stddef.h>
#include <stdint.h>

#define AES128CTR_KEYBYTES 16
#if defined(__AVX2_MODE__) || defined(__NEON_MODE__)
#define AES128CTR_BLOCKBYTES 64
#else
#define AES128CTR_BLOCKBYTES 16
#endif

// We support until 2^64 blocks
typedef struct {
#ifndef __AVX2_MODE__
  uint8_t rkeys[AES128CTR_KEYBYTES * 11];
#else
  __m128i rkeys[11];
#endif
  uint64_t n;
} aes128ctr_ctx;

// We support until 2^64 blocks
void aes128ctr_init(aes128ctr_ctx *state, const uint8_t key[AES128CTR_KEYBYTES],
                    uint64_t nonce);

// We support until 2^64 blocks
void aes128ctr_squeeze(uint8_t *out, size_t outlen, aes128ctr_ctx *state,
                       uint64_t domain_separator);

// We support until 2^64 blocks
void aes128ctr_squeeze_omp(uint8_t *out, size_t outlen, aes128ctr_ctx *state,
                           uint64_t domain_separator);

void aes128ctr_state_update(aes128ctr_ctx *state, uint64_t blocks);

#endif  // __AES128CTR_H__
