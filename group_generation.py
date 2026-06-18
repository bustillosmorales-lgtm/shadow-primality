#!/usr/bin/env python3
"""
group_generation.py

Empirical (bounded) claims about how the eight shadow bases sit inside (Z/pZ)*:

(A) The eight shadow bases sigma_r(p) generate the full group (Z/pZ)* for every
    prime 5 <= p below the first deficient prime; the paper states "1,227 primes,
    all subgroup index 1; first failure p = 6581".  We recompute the first
    deficient prime (smallest p with I_p > 1) and the exact count of primes
    p in [5, that) -- so the code reports the true (count, first-failure) pair.

(B) For p = 5 and p = 7 the shadow-base set {-r/24 mod p : r in PHI24} covers
    ALL of Z/pZ (5 distinct values mod 5, 7 mod 7); hence no composite with a
    factor 5 or 7 can have all eight bases as liars.

(C) The seven quadratic conditions of the Legendre-alignment obstruction realize
    all 2^7 = 128 sign patterns.  Paper: "verified among 9,583 primes to 10^5".
    We scan primes p <= 10^5 (p > 23), record the vector
    ( (r/p) )_{r in {5,7,11,13,17,19,23}}, and count distinct patterns (want 128)
    and the number of primes scanned.

Pure standard library.  Runtime: a few seconds.
"""
from math import gcd, isqrt
from shadow_core import is_prime, deficiency_index, PHI24


def primes_upto(B):
    sv = bytearray([1]) * (B + 1)
    sv[0] = sv[1] = 0
    for i in range(2, isqrt(B) + 1):
        if sv[i]:
            sv[i * i::i] = bytearray(len(sv[i * i::i]))
    return [i for i in range(2, B + 1) if sv[i]]


def legendre(a, p):
    a %= p
    if a == 0:
        return 0
    return 1 if pow(a, (p - 1) // 2, p) == 1 else -1


# ---------------------------------------------------------------- (A)
print("=== (A) shadow bases generate (Z/pZ)* below the first deficient prime ===")
first_def = None
count_full = 0
p = 5
while True:
    if is_prime(p) and gcd(p, 24) == 1:
        if deficiency_index(p) == 1:
            count_full += 1
        else:
            first_def = p
            break
    p += 1
print(f"  first deficient prime (smallest p with I_p > 1) : {first_def}")
print(f"  primes 5 <= p < {first_def} with full generation : {count_full}")
print(f"  paper states                                     : first failure 6581, 1227 primes")
# exact prime count below 6581, for comparison with the stated 1227
below_6581 = sum(1 for q in primes_upto(6581 - 1) if q >= 5)
print(f"  (sanity) number of primes 5 <= p < 6581          : {below_6581}")

# ---------------------------------------------------------------- (B)
print("\n=== (B) shadow bases cover all residues mod 5 and mod 7 ===")
for p in (5, 7):
    inv24 = pow(24, -1, p)
    vals = sorted(set((-r * inv24) % p for r in PHI24))
    print(f"  p={p}: shadow-base residues = {vals}  covers all of Z/{p}Z: {len(vals) == p}")

# ---------------------------------------------------------------- (C)
print("\n=== (C) all 128 Legendre patterns realized among primes <= 10^5 ===")
RAYS7 = (5, 7, 11, 13, 17, 19, 23)
patterns = set()
scanned = 0
for p in primes_upto(10**5):
    if p <= 23:
        continue
    scanned += 1
    patterns.add(tuple(legendre(r, p) for r in RAYS7))
print(f"  primes scanned (23 < p <= 10^5) : {scanned}")
print(f"  distinct Legendre patterns      : {len(patterns)}  (want 128 = 2^7)")
print(f"  all 128 realized                : {len(patterns) == 128}")
