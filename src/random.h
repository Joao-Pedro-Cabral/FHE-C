
#ifndef __RANDOM_H__
#define __RANDOM_H__

#include <stddef.h>
#include <stdint.h>

// n random bytes
// Since this function is not that relevant to the performance
// For now, I just do a memset
int randombytes(uint8_t* buf, size_t n);

#endif
