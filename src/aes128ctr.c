/* Based heavily on public-domain code by Romain Dolbeau
 * Different handling of nonce+counter than original version using
 * separated 64-bit nonce and internal 64-bit counter, starting from zero
 * Public Domain */

#include "aes128ctr.h"

#include <stdio.h>
#include <string.h>

#ifdef __NEON_MODE__
static void aesni_encrypt(uint8_t out[AES128CTR_BLOCKBYTES], uint64_t n,
                          const uint8_t rkeys[AES128CTR_KEYBYTES * 11],
                          uint64_t domain_separator) {
  uint8x16_t f, f0, f1, f2, f3;

  /* Increase counter in 4 consecutive blocks */
  f0 = {n + 0, domain_separator};
  f1 = {n + 1, domain_separator};
  f2 = {n + 2, domain_separator};
  f3 = {n + 3, domain_separator};

  /* Write counter for next iteration, increased by 4 */
  // vst1_u64((uint64_t *)n,
  //          vadd_u64(vget_low_u64((uint64x2_t)f), (uint64x1_t){4}));

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
#elif __AVX2_MODE__
static void aesni_encrypt(uint8_t out[AES128CTR_BLOCKBYTES], uint64_t n,
                          const __m128i rkeys[11], uint64_t domain_separator) {
  __m128i f, f0, f1, f2, f3;

  /* Increase counter in 4 consecutive blocks */
  f0 = _mm_set_epi64x(domain_separator, (int64_t)n);
  f1 = _mm_set_epi64x(domain_separator, (int64_t)n + 1);
  f2 = _mm_set_epi64x(domain_separator, (int64_t)n + 2);
  f3 = _mm_set_epi64x(domain_separator, (int64_t)n + 3);

  /* Write counter for next iteration, increased by 4 */
  // _mm_store_si128(n, _mm_add_epi64(f, _mm_set_epi64x(0, 4)));

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
#else
static const uint8_t sbox[256] = {
    0x63, 0x7c, 0x77, 0x7b, 0xf2, 0x6b, 0x6f, 0xc5, 0x30, 0x01, 0x67, 0x2b,
    0xfe, 0xd7, 0xab, 0x76, 0xca, 0x82, 0xc9, 0x7d, 0xfa, 0x59, 0x47, 0xf0,
    0xad, 0xd4, 0xa2, 0xaf, 0x9c, 0xa4, 0x72, 0xc0, 0xb7, 0xfd, 0x93, 0x26,
    0x36, 0x3f, 0xf7, 0xcc, 0x34, 0xa5, 0xe5, 0xf1, 0x71, 0xd8, 0x31, 0x15,
    0x04, 0xc7, 0x23, 0xc3, 0x18, 0x96, 0x05, 0x9a, 0x07, 0x12, 0x80, 0xe2,
    0xeb, 0x27, 0xb2, 0x75, 0x09, 0x83, 0x2c, 0x1a, 0x1b, 0x6e, 0x5a, 0xa0,
    0x52, 0x3b, 0xd6, 0xb3, 0x29, 0xe3, 0x2f, 0x84, 0x53, 0xd1, 0x00, 0xed,
    0x20, 0xfc, 0xb1, 0x5b, 0x6a, 0xcb, 0xbe, 0x39, 0x4a, 0x4c, 0x58, 0xcf,
    0xd0, 0xef, 0xaa, 0xfb, 0x43, 0x4d, 0x33, 0x85, 0x45, 0xf9, 0x02, 0x7f,
    0x50, 0x3c, 0x9f, 0xa8, 0x51, 0xa3, 0x40, 0x8f, 0x92, 0x9d, 0x38, 0xf5,
    0xbc, 0xb6, 0xda, 0x21, 0x10, 0xff, 0xf3, 0xd2, 0xcd, 0x0c, 0x13, 0xec,
    0x5f, 0x97, 0x44, 0x17, 0xc4, 0xa7, 0x7e, 0x3d, 0x64, 0x5d, 0x19, 0x73,
    0x60, 0x81, 0x4f, 0xdc, 0x22, 0x2a, 0x90, 0x88, 0x46, 0xee, 0xb8, 0x14,
    0xde, 0x5e, 0x0b, 0xdb, 0xe0, 0x32, 0x3a, 0x0a, 0x49, 0x06, 0x24, 0x5c,
    0xc2, 0xd3, 0xac, 0x62, 0x91, 0x95, 0xe4, 0x79, 0xe7, 0xc8, 0x37, 0x6d,
    0x8d, 0xd5, 0x4e, 0xa9, 0x6c, 0x56, 0xf4, 0xea, 0x65, 0x7a, 0xae, 0x08,
    0xba, 0x78, 0x25, 0x2e, 0x1c, 0xa6, 0xb4, 0xc6, 0xe8, 0xdd, 0x74, 0x1f,
    0x4b, 0xbd, 0x8b, 0x8a, 0x70, 0x3e, 0xb5, 0x66, 0x48, 0x03, 0xf6, 0x0e,
    0x61, 0x35, 0x57, 0xb9, 0x86, 0xc1, 0x1d, 0x9e, 0xe1, 0xf8, 0x98, 0x11,
    0x69, 0xd9, 0x8e, 0x94, 0x9b, 0x1e, 0x87, 0xe9, 0xce, 0x55, 0x28, 0xdf,
    0x8c, 0xa1, 0x89, 0x0d, 0xbf, 0xe6, 0x42, 0x68, 0x41, 0x99, 0x2d, 0x0f,
    0xb0, 0x54, 0xbb, 0x16};

static const uint8_t rcon[10] = {0x01, 0x02, 0x04, 0x08, 0x10,
                                 0x20, 0x40, 0x80, 0x1b, 0x36};

static uint8_t xtime(uint8_t x) {
  return (uint8_t)((x << 1) ^ ((x >> 7) * 0x1b));
}

static void sub_bytes(uint8_t s[16]) {
  for (int i = 0; i < 16; ++i) s[i] = sbox[s[i]];
}

static void shift_rows(uint8_t s[16]) {
  uint8_t t[16];

  /*
   * AES state is column-major:
   *
   *  0  4  8 12
   *  1  5  9 13
   *  2  6 10 14
   *  3  7 11 15
   */

  t[0] = s[0];
  t[4] = s[4];
  t[8] = s[8];
  t[12] = s[12];

  t[1] = s[5];
  t[5] = s[9];
  t[9] = s[13];
  t[13] = s[1];

  t[2] = s[10];
  t[6] = s[14];
  t[10] = s[2];
  t[14] = s[6];

  t[3] = s[15];
  t[7] = s[3];
  t[11] = s[7];
  t[15] = s[11];

  memcpy(s, t, 16);
}

static void mix_columns(uint8_t s[16]) {
  for (int c = 0; c < 4; ++c) {
    uint8_t *p = &s[4 * c];

    uint8_t a = p[0];
    uint8_t b = p[1];
    uint8_t c0 = p[2];
    uint8_t d = p[3];

    uint8_t x = a ^ b ^ c0 ^ d;

    p[0] ^= x ^ xtime(a ^ b);
    p[1] ^= x ^ xtime(b ^ c0);
    p[2] ^= x ^ xtime(c0 ^ d);
    p[3] ^= x ^ xtime(d ^ a);
  }
}

static void add_round_key(uint8_t s[16], const uint8_t rk[16]) {
  for (int i = 0; i < 16; ++i) s[i] ^= rk[i];
}

static void add_nounce(uint8_t s[16], uint64_t n, uint64_t domain_separator) {
  uint64_t tmp = n;
  for (int i = 0; i < 16; ++i) {
    s[i] = (uint8_t)tmp;
    tmp >>= 8;
    if (i == 8) tmp = domain_separator;
  }
}

static void aesni_encrypt(uint8_t out[AES128CTR_BLOCKBYTES], uint64_t n,
                          const uint8_t rkeys[AES128CTR_KEYBYTES * 11],
                          uint64_t domain_separator) {
  uint8_t s[16];

  add_nounce(s, n, domain_separator);

  add_round_key(s, &rkeys[0]);

  for (int round = 1; round < 10; ++round) {
    sub_bytes(s);
    shift_rows(s);
    mix_columns(s);
    add_round_key(s, &rkeys[AES128CTR_KEYBYTES * round]);
  }

  sub_bytes(s);
  shift_rows(s);
  add_round_key(s, &rkeys[AES128CTR_KEYBYTES * 10]);

  memcpy(out, s, 16);
}
#endif

void aes128ctr_init(aes128ctr_ctx *state, const uint8_t key[AES128CTR_KEYBYTES],
                    uint64_t nonce) {
  state->n = nonce;
#if defined(__NEON_MODE__) || defined(__AVX2_MODE__)
  int idx = 0;
#endif
#ifdef __NEON_MODE__
  uint8x16_t key0, zeros, s8, s8_2, imm;

  key0 = vld1q_u8(key);
  uint8x16_t idxs = {9, 6, 3, 12, 9, 6, 3, 12, 9, 6, 3, 12, 9, 6, 3, 12};

  s8 = key0;
  zeros = vdupq_n_u8(0);

#define BLOCK1(IMM)                                      \
  vst1q_u8(&state->rkeys[AES128CTR_KEYBYTES * idx], s8); \
  idx++;                                                 \
  imm = vsetq_lane_u8(IMM, zeros, 9);                    \
  s8_2 = veorq_u8(vaeseq_u8(s8, zeros), imm);            \
  s8_2 = vqtbl1q_u8(s8_2, idxs);                         \
  s8 = veorq_u8(s8, vextq_u8(zeros, s8, 12));            \
  s8 = veorq_u8(s8, vextq_u8(zeros, s8, 8));             \
  s8 = veorq_u8(s8, s8_2);

#elif __AVX2_MODE__

  __m128i key0, temp0, temp1, temp4;

  key0 = _mm_loadu_si128((__m128i *)(key));

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

#if defined(__NEON_MODE__) || defined(__AVX2_MODE__)

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

#endif

#ifdef __NEON_MODE__
  vst1q_u8(&state->rkeys[AES128CTR_KEYBYTES * idx], s8);
#elif __AVX2_MODE__
  state->rkeys[idx++] = temp0;
#else
  /* Round 0 */
  for (int i = 0; i < 16; i++) state->rkeys[i] = key[i];

  for (int round = 1; round <= 10; round++) {
    const uint8_t *prev = &state->rkeys[AES128CTR_KEYBYTES * (round - 1)];
    uint8_t *rk = &state->rkeys[AES128CTR_KEYBYTES * round];

    /* First word = SubWord(RotWord(previous last word)) ^ Rcon */
    rk[0] = prev[0] ^ sbox[prev[13]] ^ rcon[round - 1];
    rk[1] = prev[1] ^ sbox[prev[14]];
    rk[2] = prev[2] ^ sbox[prev[15]];
    rk[3] = prev[3] ^ sbox[prev[12]];

    /* Remaining words */
    for (int i = 4; i < 16; i++) rk[i] = prev[i] ^ rk[i - 4];
  }
#endif
}

void aes128ctr_squeeze(uint8_t *out, size_t outlen, aes128ctr_ctx *state,
                       size_t domain_separator) {
  uint8_t buf[AES128CTR_BLOCKBYTES];
  uint64_t n = state->n;

  while (outlen >= AES128CTR_BLOCKBYTES) {
#ifdef __NEON_MODE__
    aesni_encrypt(out, n, state->rkeys, domain_separator);
#elif __AVX2_MODE__
    aesni_encrypt(out, n, state->rkeys, domain_separator);
#else
    aesni_encrypt(out, n, state->rkeys, domain_separator);
#endif
    n += 4;
    outlen -= AES128CTR_BLOCKBYTES;
    out += AES128CTR_BLOCKBYTES;
  }

  if (outlen) {
#ifdef __NEON_MODE__
    aesni_encrypt(buf, n, state->rkeys, domain_separator);
#elif __AVX2_MODE__
    aesni_encrypt(buf, n, state->rkeys, domain_separator);
#else
    aesni_encrypt(buf, n, state->rkeys, domain_separator);
#endif
    n += 4;
    for (size_t i = 0; i < outlen; i++) out[i] = buf[i];
  }
}

void aes128ctr_squeeze_omp(uint8_t *out, size_t outlen, aes128ctr_ctx *state,
                           size_t domain_separator) {
  uint8_t buf[AES128CTR_BLOCKBYTES];
  uint64_t blocks = outlen / AES128CTR_BLOCKBYTES;
  uint64_t remaining = outlen - blocks * AES128CTR_BLOCKBYTES;
  uint64_t n = state->n;

#pragma omp for nowait
  for (uint64_t i = 0; i < blocks; i++) {
#ifdef __NEON_MODE__
    aesni_encrypt(&out[AES128CTR_BLOCKBYTES * i], n + 4 * i, state->rkeys,
                  domain_separator);
#elif __AVX2_MODE__
    aesni_encrypt(&out[AES128CTR_BLOCKBYTES * i], n + 4 * i, state->rkeys,
                  domain_separator);
#else
    aesni_encrypt(&out[AES128CTR_BLOCKBYTES * i], n + 4 * i, state->rkeys,
                  domain_separator);
#endif
  }

#pragma omp single
  {
    // state->n = state->n + 4 * (blocks + remaining);
    n = n + 4 * (blocks + remaining);
    if (remaining > 0) {
#ifdef __NEON_MODE__
      aesni_encrypt(buf, n, state->rkeys, domain_separator);
#elif __AVX2_MODE__
      aesni_encrypt(buf, n, state->rkeys, domain_separator);
#else
      aesni_encrypt(buf, n, state->rkeys, domain_separator);
#endif
      for (size_t i = 0; i < remaining; i++)
        out[blocks * AES128CTR_BLOCKBYTES + i] = buf[i];
    }
  }
}

void aes128ctr_state_update(aes128ctr_ctx *state, uint64_t blocks) {
  state->n += blocks;
}
