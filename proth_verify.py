#!/usr/bin/env python3
"""
Verify the 10,006-digit certified Proth prime claimed in the paper:
    N = 526623 * 2^33220 + 1.

Proth's theorem: for N = k*2^m + 1 with k odd and k < 2^m, N is prime
iff there exists an integer a with  a^((N-1)/2) = -1 (mod N).
Here k = 526623 (odd), m = 33220, and k < 2^m, so the criterion applies.

A single modular exponentiation on N certifies primality. We also report a
small set of bases for which the Proth witness holds.
"""
import time, sys
sys.set_int_max_str_digits(20000)

k, m = 526623, 33220
N = k * (1 << m) + 1
half = (N - 1) // 2
digits = len(str(N))

print(f"N = 526623 * 2^33220 + 1")
print(f"  digits          : {digits}")
print(f"  k odd           : {k % 2 == 1}")
print(f"  k < 2^m         : {k < (1 << m)}")
print()

t0 = time.time()
witnesses = []
for a in (3, 5, 7, 11, 13):
    x = pow(a, half, N)
    is_w = (x == N - 1)              # a^((N-1)/2) = -1 (mod N)
    print(f"  base a={a:>2}: a^((N-1)/2) mod N == N-1 ? {is_w}   ({time.time()-t0:.1f}s elapsed)")
    if is_w:
        witnesses.append(a)
        break   # one Proth witness already certifies primality

if witnesses:
    print(f"\nRESULT: N is PRIME (Proth witness a={witnesses[0]}). Certified in {time.time()-t0:.1f}s.")
else:
    print("\nRESULT: no Proth witness among tested bases (N composite or try more bases).")
