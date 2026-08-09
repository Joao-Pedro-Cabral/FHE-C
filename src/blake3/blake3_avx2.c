
#include "blake3_impl.h"
#include <immintrin.h>
#include <stdio.h>

static inline __m256i _mm256_rotr32(__m256i x, const int n)
{
  return _mm256_or_si256(_mm256_srli_epi32(x, n), _mm256_slli_epi32(x, 32 - n));
}

INLINE void g(__m256i *stateA, __m256i *stateB, __m256i *stateC, __m256i *stateD,
              uint32_t x, uint32_t y) {
  __m256i vx = _mm256_set1_epi32(x);
  __m256i vy = _mm256_set1_epi32(y);
  *stateA = _mm256_add_epi32(*stateA, _mm256_add_epi32(*stateB, vx));
  *stateD = _mm256_rotr32(_mm256_xor_si256(*stateA, *stateD), 16);
  *stateC = _mm256_add_epi32(*stateC, *stateD);
  *stateB = _mm256_rotr32(_mm256_xor_si256(*stateB, *stateC), 12);

  *stateA = _mm256_add_epi32(*stateA, _mm256_add_epi32(*stateB, vy));
  *stateD = _mm256_rotr32(_mm256_xor_si256(*stateA, *stateD), 8);
  *stateC = _mm256_add_epi32(*stateC, *stateD);
  *stateB = _mm256_rotr32(_mm256_xor_si256(*stateB, *stateC), 7);
}

INLINE void init_state(__m256i *state0,  __m256i *state1,  __m256i *state2,  __m256i *state3,
                       __m256i *state4,  __m256i *state5,  __m256i *state6,  __m256i *state7,
                       __m256i *state8,  __m256i *state9,  __m256i *state10, __m256i *state11,
                       __m256i *state12, __m256i *state13, __m256i *state14, __m256i *state15,
                       const uint32_t cv[8], uint8_t block_len, uint64_t counter, uint8_t flags) {
  static const uint32_t idx[8] __attribute__((aligned(32))) = {
      0, 1, 2, 3, 4, 5, 6, 7
  };
  __m256i inc = _mm256_load_si256((const __m256i *)idx);

  *state0 = _mm256_set1_epi32(cv[0]);
  *state1 = _mm256_set1_epi32(cv[1]);
  *state2 = _mm256_set1_epi32(cv[2]);
  *state3 = _mm256_set1_epi32(cv[3]);
  *state4 = _mm256_set1_epi32(cv[4]);
  *state5 = _mm256_set1_epi32(cv[5]);
  *state6 = _mm256_set1_epi32(cv[6]);
  *state7 = _mm256_set1_epi32(cv[7]);
  *state8 = _mm256_set1_epi32(IV[0]);
  *state9 = _mm256_set1_epi32(IV[1]);
  *state10 = _mm256_set1_epi32(IV[2]);
  *state11 = _mm256_set1_epi32(IV[3]);
  *state12 = _mm256_add_epi32(_mm256_set1_epi32(counter_low(counter)), inc);
  *state13 = _mm256_set1_epi32(counter_high(counter));
  *state14 = _mm256_set1_epi32((uint32_t)block_len);
  *state15 = _mm256_set1_epi32((uint32_t)flags);
}

INLINE void round_fn(__m256i *state0,  __m256i *state1,  __m256i *state2,  __m256i *state3,
                     __m256i *state4,  __m256i *state5,  __m256i *state6,  __m256i *state7,
                     __m256i *state8,  __m256i *state9,  __m256i *state10, __m256i *state11,
                     __m256i *state12, __m256i *state13, __m256i *state14, __m256i *state15,
                     const uint32_t *msg, size_t round) {
  // Select the message schedule based on the round.
  const uint8_t *schedule = MSG_SCHEDULE[round];

  // Mix the columns.
  g(state0, state4, state8,  state12, msg[schedule[0]], msg[schedule[1]]);
  g(state1, state5, state9,  state13, msg[schedule[2]], msg[schedule[3]]);
  g(state2, state6, state10, state14, msg[schedule[4]], msg[schedule[5]]);
  g(state3, state7, state11, state15, msg[schedule[6]], msg[schedule[7]]);

  // Mix the rows.
  g(state0, state5, state10, state15, msg[schedule[8]],  msg[schedule[9]]);
  g(state1, state6, state11, state12, msg[schedule[10]], msg[schedule[11]]);
  g(state2, state7, state8,  state13, msg[schedule[12]], msg[schedule[13]]);
  g(state3, state4, state9,  state14, msg[schedule[14]], msg[schedule[15]]);
}

// If we transpose the first(second) group of 8 registers, we've the even(odd) blocks of the output
INLINE void transpose_store(__m256i *state0,  __m256i *state1,  __m256i *state2,  __m256i *state3,
                            __m256i *state4,  __m256i *state5,  __m256i *state6,  __m256i *state7,
                            uint8_t* out) {
  __m256i tmp;
  // First layer
  tmp = *state0;
  *state0 = _mm256_unpacklo_epi32(tmp, *state1);
  *state1 = _mm256_unpackhi_epi32(tmp, *state1);
  tmp = *state2;
  *state2 = _mm256_unpacklo_epi32(tmp, *state3);
  *state3 = _mm256_unpackhi_epi32(tmp, *state3);
  tmp = *state4;
  *state4 = _mm256_unpacklo_epi32(tmp, *state5);
  *state5 = _mm256_unpackhi_epi32(tmp, *state5);
  tmp = *state6;
  *state6 = _mm256_unpacklo_epi32(tmp, *state7);
  *state7 = _mm256_unpackhi_epi32(tmp, *state7);
  // Second Layer
  tmp = *state0;
  *state0 = _mm256_unpacklo_epi64(tmp, *state2);
  *state2 = _mm256_unpackhi_epi64(tmp, *state2);
  tmp = *state1;
  *state1 = _mm256_unpacklo_epi64(tmp, *state3);
  *state3 = _mm256_unpackhi_epi64(tmp, *state3);
  tmp = *state4;
  *state4 = _mm256_unpacklo_epi64(tmp, *state6);
  *state6 = _mm256_unpackhi_epi64(tmp, *state6);
  tmp = *state5;
  *state5 = _mm256_unpacklo_epi64(tmp, *state7);
  *state7 = _mm256_unpackhi_epi64(tmp, *state7);
  // Third Layer
  tmp = *state0;
  *state0 = _mm256_permute2x128_si256(tmp, *state4, 0x20);
  *state4 = _mm256_permute2x128_si256(tmp, *state4, 0x31);
  tmp = *state1;
  *state1 = _mm256_permute2x128_si256(tmp, *state5, 0x20);
  *state5 = _mm256_permute2x128_si256(tmp, *state5, 0x31);
  tmp = *state2;
  *state2 = _mm256_permute2x128_si256(tmp, *state6, 0x20);
  *state6 = _mm256_permute2x128_si256(tmp, *state6, 0x31);
  tmp = *state3;
  *state3 = _mm256_permute2x128_si256(tmp, *state7, 0x20);
  *state7 = _mm256_permute2x128_si256(tmp, *state7, 0x31);
  // Store
  _mm256_storeu_si256((__m256i*)out, *state0);
  _mm256_storeu_si256((__m256i*)(out + 64), *state2);
  _mm256_storeu_si256((__m256i*)(out + 128), *state1);
  _mm256_storeu_si256((__m256i*)(out + 192), *state3);
  _mm256_storeu_si256((__m256i*)(out + 256), *state4);
  _mm256_storeu_si256((__m256i*)(out + 320), *state6);
  _mm256_storeu_si256((__m256i*)(out + 384), *state5);
  _mm256_storeu_si256((__m256i*)(out + 448), *state7);
}

INLINE void finalize_state(__m256i *state0,  __m256i *state1,  __m256i *state2,  __m256i *state3,
                           __m256i *state4,  __m256i *state5,  __m256i *state6,  __m256i *state7,
                           __m256i *state8,  __m256i *state9,  __m256i *state10, __m256i *state11,
                           __m256i *state12, __m256i *state13, __m256i *state14, __m256i *state15,
                           const uint32_t cv[8]) {
  *state0 = _mm256_xor_si256(*state0, *state8);
  *state1 = _mm256_xor_si256(*state1, *state9);
  *state2 = _mm256_xor_si256(*state2, *state10);
  *state3 = _mm256_xor_si256(*state3, *state11);
  *state4 = _mm256_xor_si256(*state4, *state12);
  *state5 = _mm256_xor_si256(*state5, *state13);
  *state6 = _mm256_xor_si256(*state6, *state14);
  *state7 = _mm256_xor_si256(*state7, *state15);
  *state8 =  _mm256_xor_si256(*state8 , _mm256_set1_epi32(cv[0]));
  *state9 =  _mm256_xor_si256(*state9 , _mm256_set1_epi32(cv[1]));
  *state10 = _mm256_xor_si256(*state10, _mm256_set1_epi32(cv[2]));
  *state11 = _mm256_xor_si256(*state11, _mm256_set1_epi32(cv[3]));
  *state12 = _mm256_xor_si256(*state12, _mm256_set1_epi32(cv[4]));
  *state13 = _mm256_xor_si256(*state13, _mm256_set1_epi32(cv[5]));
  *state14 = _mm256_xor_si256(*state14, _mm256_set1_epi32(cv[6]));
  *state15 = _mm256_xor_si256(*state15, _mm256_set1_epi32(cv[7]));
}

void blake3_xof_many_avx2(const uint32_t cv[8],
                          const uint8_t block[BLAKE3_BLOCK_LEN],
                          uint8_t block_len, uint64_t counter, uint8_t flags,
                          uint8_t* out, size_t outblocks) {
  size_t outblocks8 = outblocks - (outblocks & 0x7);
  uint32_t block_words[16];
  __m256i state0, state1, state2, state3, state4, state5, state6, state7,
          state8, state9, state10, state11, state12, state13, state14, state15;
  if (outblocks8 > 0) {
    block_words[0] = load32(block + 4 * 0);
    block_words[1] = load32(block + 4 * 1);
    block_words[2] = load32(block + 4 * 2);
    block_words[3] = load32(block + 4 * 3);
    block_words[4] = load32(block + 4 * 4);
    block_words[5] = load32(block + 4 * 5);
    block_words[6] = load32(block + 4 * 6);
    block_words[7] = load32(block + 4 * 7);
    block_words[8] = load32(block + 4 * 8);
    block_words[9] = load32(block + 4 * 9);
    block_words[10] = load32(block + 4 * 10);
    block_words[11] = load32(block + 4 * 11);
    block_words[12] = load32(block + 4 * 12);
    block_words[13] = load32(block + 4 * 13);
    block_words[14] = load32(block + 4 * 14);
    block_words[15] = load32(block + 4 * 15);
  }
  // We process 256/32 = 8 blocks per time
  for (size_t i = 0; i < outblocks8; i += 8) {
    init_state(&state0, &state1, &state2, &state3, &state4, &state5, &state6, &state7,
               &state8, &state9, &state10, &state11, &state12, &state13, &state14, &state15,
              cv, block_len, counter, flags);
    round_fn(&state0, &state1, &state2, &state3, &state4, &state5, &state6, &state7,
             &state8, &state9, &state10, &state11, &state12, &state13, &state14, &state15,
             block_words, 0);
    round_fn(&state0, &state1, &state2, &state3, &state4, &state5, &state6, &state7,
             &state8, &state9, &state10, &state11, &state12, &state13, &state14, &state15,
             block_words, 1);
    round_fn(&state0, &state1, &state2, &state3, &state4, &state5, &state6, &state7,
             &state8, &state9, &state10, &state11, &state12, &state13, &state14, &state15,
             block_words, 2);
    round_fn(&state0, &state1, &state2, &state3, &state4, &state5, &state6, &state7,
             &state8, &state9, &state10, &state11, &state12, &state13, &state14, &state15,
             block_words, 3);
    round_fn(&state0, &state1, &state2, &state3, &state4, &state5, &state6, &state7,
             &state8, &state9, &state10, &state11, &state12, &state13, &state14, &state15,
             block_words, 4);
    round_fn(&state0, &state1, &state2, &state3, &state4, &state5, &state6, &state7,
             &state8, &state9, &state10, &state11, &state12, &state13, &state14, &state15,
             block_words, 5);
    round_fn(&state0, &state1, &state2, &state3, &state4, &state5, &state6, &state7,
             &state8, &state9, &state10, &state11, &state12, &state13, &state14, &state15,
             block_words, 6);
    finalize_state(&state0, &state1, &state2, &state3, &state4, &state5, &state6, &state7,
                   &state8, &state9, &state10, &state11, &state12, &state13, &state14, &state15,
                   cv);
    transpose_store(&state0, &state1, &state2, &state3, &state4, &state5, &state6, &state7, out + 64*i);
    transpose_store(&state8, &state9, &state10, &state11, &state12, &state13, &state14, &state15, out + 64*i + 32);
    counter += 8;
  }
  for (size_t i = outblocks8; i < outblocks; ++i) {
    blake3_compress_xof(cv, block, block_len, counter + i, flags, out + 64*i);
  }
}
