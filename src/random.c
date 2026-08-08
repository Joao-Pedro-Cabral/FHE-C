
#include "random.h"

#include "string.h"

static uint8_t randomness = 0;

void randombytes(uint8_t* buf, uint32_t n) {
  memset(buf, randomness, n);
  randomness++;
}
