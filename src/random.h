
#ifndef __RANDOM_H__
#define __RANDOM_H__

#include <stdint.h>

// n random bytes
// Since this function is not that relevant to the performance
// For now, I just do a memset
void randombytes(uint8_t* buf, uint32_t n);

#endif
