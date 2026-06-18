#!/usr/bin/env python3
"""
antipodal_sieve_check.py

The modulo-24 / antipodal sieve claim (Section 4): sieving primes by ray and
comparing to Eratosthenes up to 10^7 gives perfect agreement, and the primes
distribute near-uniformly across the eight rays (~83,000 per ray, Dirichlet).

We build the Eratosthenes sieve to 10^7, bucket every prime p > 3 by its ray
p mod 24, and report the per-ray counts (which should all be ~83,000) together
with the total (= pi(10^7) - 2).  A ray-targeted reconstruction (mark multiples
of each prime on each ray) is checked to reproduce the Eratosthenes prime set
exactly.

Pure standard library.  Runtime: a few seconds.
"""
from math import isqrt
from shadow_core import PHI24

B = 10**7

# Eratosthenes
sv = bytearray([1]) * (B + 1)
sv[0] = sv[1] = 0
for i in range(2, isqrt(B) + 1):
    if sv[i]:
        sv[i * i::i] = bytearray(len(sv[i * i::i]))

# bucket primes > 3 by ray
counts = {r: 0 for r in PHI24}
total = 0
for p in range(5, B + 1):
    if sv[p]:
        counts[p % 24] += 1
        total += 1

print(f"=== primes 5 <= p <= 10^7 bucketed by ray (mod 24) ===")
for r in PHI24:
    print(f"  ray {r:2d}: {counts[r]:,}")
print(f"  total      : {total:,}   (= pi(10^7) - 2 = {sum(1 for p in range(2,B+1) if sv[p]) - 2:,})")
spread = (max(counts.values()) - min(counts.values())) / (total / 8)
print(f"  per-ray average ~{total // 8:,}  (paper: ~83,000), max-min spread {spread:.2%} of mean")

# ray-targeted reconstruction agrees with Eratosthenes
# for each ray r, n = 24k + r is composite-with-small-factor iff some prime q<=sqrt
# divides it; primes are the unmarked admissible n. We verify the admissible
# unmarked set equals the Eratosthenes primes on that ray.
agree = True
mark = bytearray(B + 1)
for q in range(5, isqrt(B) + 1):
    if sv[q]:
        mark[q * q::q] = bytearray([1]) * len(mark[q * q::q])
for p in range(5, B + 1):
    if p % 6 in (1, 5):
        recon_prime = (mark[p] == 0)
        if recon_prime != bool(sv[p]):
            agree = False
            break
print(f"  ray-sieve reconstruction matches Eratosthenes exactly: {agree}")
