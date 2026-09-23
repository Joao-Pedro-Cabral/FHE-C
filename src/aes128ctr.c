/* Based heavily on public-domain code by Romain Dolbeau
 * Different handling of nonce+counter than original version using
 * separated 64-bit nonce and internal 64-bit counter, starting from zero
 * Public Domain */

#include "aes128ctr.h"

#include <stdio.h>

#ifdef __aarch64__
static void aesni_encrypt4(uint8_t out[AES128CTR_BLOCKBYTES],
                           uint8_t n[AES128CTR_KEYBYTES],
                           const uint8_t rkeys[AES128CTR_KEYBYTES * 11]) {
  uint8x16_t f, f0, f1, f2, f3;

  /* Load current counter value */
  f = vld1q_u8(n);

  /* Increase counter in 4 consecutive blocks */
  f0 = vrev64q_u8((uint8x16_t)vaddq_u64((uint64x2_t)f, (uint64x2_t){0, 0}));
  f1 = vrev64q_u8((uint8x16_t)vaddq_u64((uint64x2_t)f, (uint64x2_t){1, 0}));
  f2 = vrev64q_u8((uint8x16_t)vaddq_u64((uint64x2_t)f, (uint64x2_t){2, 0}));
  f3 = vrev64q_u8((uint8x16_t)vaddq_u64((uint64x2_t)f, (uint64x2_t){3, 0}));

  /* Write counter for next iteration, increased by 4 */
  vst1_u64((uint64_t *)n,
           vadd_u64(vget_low_u64((uint64x2_t)f), (uint64x1_t){3}));

  /* Actual AES encryption, 4x interleaved */
  for (int i = 0; i < 9; i++) {
    f = vld1q_u8(&rkeys[AES128CTR_KEYBYTES * i]);
    f0 = vaesmcq_u8(vaeseq_u8(f0, f));
    f1 = vaesmcq_u8(vaeseq_u8(f1, f));
    f2 = vaesmcq_u8(vaeseq_u8(f2, f));
    f3 = vaesmcq_u8(vaeseq_u8(f3, f));
  }

  f = vld1q_u8(&rkeys[AES128CTR_KEYBYTES * 9]);
  f0 = vaeseq_u8(f0, f);
  f1 = vaeseq_u8(f1, f);
  f2 = vaeseq_u8(f2, f);
  f3 = vaeseq_u8(f3, f);

  f = vld1q_u8(&rkeys[AES128CTR_KEYBYTES * 10]);
  f0 = veorq_u8(f0, f);
  f1 = veorq_u8(f1, f);
  f2 = veorq_u8(f2, f);
  f3 = veorq_u8(f3, f);

  /* Write results */
  vst1q_u8((out + 0), f0);
  vst1q_u8((out + AES128CTR_KEYBYTES), f1);
  vst1q_u8((out + 2 * AES128CTR_KEYBYTES), f2);
  vst1q_u8((out + 3 * AES128CTR_KEYBYTES), f3);
}
#else
static void aesni_encrypt4(uint8_t out[AES128CTR_BLOCKBYTES], __m128i *n,
                           const __m128i rkeys[11]) {
  __m128i f, f0, f1, f2, f3;
  const __m128i idx =
      _mm_set_epi8(8, 9, 10, 11, 12, 13, 14, 15, 7, 6, 5, 4, 3, 2, 1, 0);

  /* Load current counter value */
  f = _mm_load_si128(n);

  /* Increase counter in 4 consecutive blocks */
  f0 = _mm_shuffle_epi8(_mm_add_epi64(f, _mm_set_epi64x(0, 0)), idx);
  f1 = _mm_shuffle_epi8(_mm_add_epi64(f, _mm_set_epi64x(0, 1)), idx);
  f2 = _mm_shuffle_epi8(_mm_add_epi64(f, _mm_set_epi64x(0, 2)), idx);
  f3 = _mm_shuffle_epi8(_mm_add_epi64(f, _mm_set_epi64x(0, 3)), idx);

  /* Write counter for next iteration, increased by 4 */
  _mm_store_si128(n, _mm_add_epi64(f, _mm_set_epi64x(0, 4)));

  /* Actual AES encryption, 4x interleaved */
  f = _mm_load_si128(&rkeys[0]);
  f0 = _mm_xor_si128(f0, f);
  f1 = _mm_xor_si128(f1, f);
  f2 = _mm_xor_si128(f2, f);
  f3 = _mm_xor_si128(f3, f);

  for (int i = 1; i < 10; i++) {
    f = _mm_load_si128(&rkeys[i]);
    f0 = _mm_aesenc_si128(f0, f);
    f1 = _mm_aesenc_si128(f1, f);
    f2 = _mm_aesenc_si128(f2, f);
    f3 = _mm_aesenc_si128(f3, f);
  }

  f = _mm_load_si128(&rkeys[10]);
  f0 = _mm_aesenclast_si128(f0, f);
  f1 = _mm_aesenclast_si128(f1, f);
  f2 = _mm_aesenclast_si128(f2, f);
  f3 = _mm_aesenclast_si128(f3, f);

  /* Write results */
  _mm_storeu_si128((__m128i *)(out + 0), f0);
  _mm_storeu_si128((__m128i *)(out + 16), f1);
  _mm_storeu_si128((__m128i *)(out + 32), f2);
  _mm_storeu_si128((__m128i *)(out + 48), f3);
}
#endif

void aes128ctr_init(aes128ctr_ctx *state, const uint8_t key[AES128CTR_KEYBYTES],
                    uint64_t nonce) {
  int idx = 0;
#ifdef __aarch64__
  uint8x16_t key0, temp0, temp1, temp2, temp3, temp4, zeros;

  uint8x16_t s8, s8_2, imm;
  uint32x4_t s32, s32_2;

  key0 = vld1q_u8(key);
  uint64x2_t vn = {nonce, 0};
  uint32x4_t shifts = {-24, -16, -8, 0};
  uint32x4_t mask8 = vdupq_n_u32(0x000000FF);
  uint32x4_t mask16 = vdupq_n_u32(0x0000FF00);
  uint32x4_t mask24 = vdupq_n_u32(0x00FF0000);
  vst1q_u64((uint64_t *)state->n, vn);

  s8 = key0;
  zeros = vdupq_n_u8(0);
  temp2 = zeros;

#define BLOCK1(IMM)                                         \
  vst1q_u8(&state->rkeys[AES128CTR_KEYBYTES * idx], s8);    \
  idx++;                                                    \
  s32 = (uint32x4_t)s8;                                     \
  s32 = vextq_u32(s32, s32, 1);                             \
  imm = vdupq_n_u8(IMM);                                    \
  s8_2 = (uint8x16_t)s32;                                   \
  s8_2 = veorq_u8(vaeseq_u8(s8_2, zeros), imm);             \
  s32_2 = (uint32x4_t)s8_2;                                 \
  s32_2 = vandq_u32(vshlq_u32(s32_2, shifts), mask8);       \
  s32 = (uint32x4_t)s8;                                     \
  s32 = veorq_u32(s32, s32_2);                              \
  s32_2 = vshlq_u32(vandq_u32(s32, mask8), 8);              \
  s32 = veorq_u32(s32, s32_2);                              \
  s32_2 = vshlq_u32(vandq_u32(s32, mask16), 8);             \
  s32 = veorq_u32(s32, s32_2);                              \
  s32_2 = vshlq_u32(vandq_u32(s32, mask24), 8);             \
  s32 = veorq_u32(s32, s32_2);                              \
  s8 = (uint8x16_t)s32;                                     \

#else

  __m128i key0, temp0, temp1, temp4;

  key0 = _mm_loadu_si128((__m128i *)(key));
  state->n = _mm_loadl_epi64((__m128i *)&nonce);

  temp0 = key0;
  temp4 = _mm_setzero_si128();

#define BLOCK1(IMM)                                                    \
  temp1 = _mm_aeskeygenassist_si128(temp0, IMM);                       \
  state->rkeys[idx++] = temp0;                                         \
  temp4 = (__m128i)_mm_shuffle_ps((__m128)temp4, (__m128)temp0, 0x10); \
  temp0 = _mm_xor_si128(temp0, temp4);                                 \
  temp4 = (__m128i)_mm_shuffle_ps((__m128)temp4, (__m128)temp0, 0x8c); \
  temp0 = _mm_xor_si128(temp0, temp4);                                 \
  temp1 = (__m128i)_mm_shuffle_ps((__m128)temp1, (__m128)temp1, 0xff); \
  temp0 = _mm_xor_si128(temp0, temp1)

#endif

  BLOCK1(0x01);

  BLOCK1(0x02);

  BLOCK1(0x04);

  BLOCK1(0x08);

  BLOCK1(0x10);

  BLOCK1(0x20);

  BLOCK1(0x40);

  BLOCK1(0x80);

  BLOCK1(0x1B);

  BLOCK1(0x36);

#ifdef __aarch64__
  vst1q_u8(&state->rkeys[AES128CTR_KEYBYTES * idx], s8);
#else
  state->rkeys[idx++] = temp0;
#endif
}

void aes128ctr_squeezeblocks(uint8_t *out, size_t nblocks,
                             aes128ctr_ctx *state) {
  for (size_t i = 0; i < nblocks; i++) {
#ifdef __aarch64__
    aesni_encrypt4(out, state->n, state->rkeys);
#else
    aesni_encrypt4(out, &state->n, state->rkeys);
#endif
    out += AES128CTR_BLOCKBYTES;
  }
}

void aes128ctr_prf(uint8_t *out, size_t outlen,
                   const uint8_t seed[AES128CTR_KEYBYTES], uint64_t nonce) {
  uint8_t buf[AES128CTR_BLOCKBYTES];
  aes128ctr_ctx state;

  aes128ctr_init(&state, seed, nonce);

  while (outlen >= AES128CTR_BLOCKBYTES) {
#ifdef __aarch64__
    aesni_encrypt4(out, state.n, state.rkeys);
#else
    aesni_encrypt4(out, &state.n, state.rkeys);
#endif
    outlen -= AES128CTR_BLOCKBYTES;
    out += AES128CTR_BLOCKBYTES;
  }

  if (outlen) {
#ifdef __aarch64__
    aesni_encrypt4(buf, state.n, state.rkeys);
#else
    aesni_encrypt4(buf, &state.n, state.rkeys);
#endif
    for (size_t i = 0; i < outlen; i++) out[i] = buf[i];
  }
}
