#!/usr/bin/env python3

import math

N = 21
SCALE = 1 << 63

real = []
imag = []
real_q63 = []
imag_q63 = []

for k in range(1, N):
    n = 1 << k
    theta = math.pi / n

    r = math.cos(theta)
    i = math.sin(theta)

    real.append(r)
    imag.append(i)

    # Saturate +1.0, since it cannot be represented in Q1.63.
    if r >= 1.0:
        rq = (1 << 63) - 1
    else:
        rq = int(round(r * SCALE))

    if i >= 1.0:
        iq = (1 << 63) - 1
    else:
        iq = int(round(i * SCALE))

    real_q63.append(rq)
    imag_q63.append(iq)


def print_double_array(name, values):
    print(f"static const double {name}[{len(values)}] = {{")
    for v in values:
        print(f"    {v:.20f},")
    print("};\n")


def print_q63_array(name, values):
    print(f"static const int64_t {name}[{len(values)}] = {{")
    for v in values:
        print(f"    INT64_C({v}),")
    print("};\n")


print("#include <stdint.h>\n")

print_double_array("twiddle_real", real)
print_double_array("twiddle_imag", imag)

print_q63_array("twiddle_real_q63", real_q63)
print_q63_array("twiddle_imag_q63", imag_q63)
