#ifndef FIPS202_H
#define FIPS202_H

#include <stddef.h>
#include <stdint.h>

#define SHAKE128_RATE 168

typedef struct {
  uint64_t ctx[25];
} shake128ctx;

/* One-stop SHAKE128 call */
void shake128(uint8_t *output, size_t outlen, const uint8_t *input,
              size_t inlen);

/* Initialize the state and absorb the provided input.
 *
 * This function does not support being called multiple times
 * with the same state.
 */
void shake128_absorb(shake128ctx *state, const uint8_t *input, size_t inlen);
/* Squeeze output out of the sponge.
 *
 * Supports being called multiple times
 */
void shake128_squeezeblocks(uint8_t *output, size_t nblocks, shake128ctx *state);

#endif