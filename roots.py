import numpy as np


def brv(x):
    """
    Returns a copy of x with its elements in bit-reversed order.
    len(x) must be a power of two.
    """
    n = len(x)

    if n == 0 or (n & (n - 1)):
        raise ValueError("Array length must be a power of two.")

    bits = n.bit_length() - 1
    y = np.empty_like(x)

    for i in range(n):
        # Reverse the lowest 'bits' bits of i
        j = int(f"{i:0{bits}b}"[::-1], 2)
        y[j] = x[i]

    return y


n = 4
k = 2

root = np.exp(1j * np.pi / 4)

x = np.array([1+4j, 3-1j, 3+1j, 1-4j], dtype=complex)

w = brv([pow(root, i) for i in range(n)])

print(x)

print(w)

for l in range(k):
    for i in range(2**(k-1-l)):
        phi = w[i]
        phi = phi.conjugate()
        for j in range(i*(2**(l+1)), i*(2**(l+1)) + (2**l)):
            t0 = x[j] + x[j + (2**l)]
            t1 = x[j] - x[j + (2**l)]
            x[j] = t0
            x[j + (2**l)] = t1*phi

print(x)
