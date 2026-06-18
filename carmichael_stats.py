#!/usr/bin/env python3
"""
carmichael_stats.py

Two bounded liar-census claims of the paper, recomputed from scratch.

PART 1 -- Carmichael numbers coprime to six up to 10^7.
  Paper: there are 102 of them; 100% have at least one shadow WITNESS
  (average 7.75 of 8 witnesses, i.e. average 0.25 liars); the best Legendre
  agreement across factors is 4 of the 7 rays, attained at
  252601 = 41 * 61 * 101.
  (coprime to 6 == coprime to 24, so all are in Phi_24.)

PART 2 -- squarefree composites with >= 3 prime factors, coprime to six, up to
  2*10^6.  Paper: 184,510 of them; ZERO have all eight strong liars; only 130
  have even one strong liar.

Method: a prime sieve + a Fermat base-2 filter isolates Carmichael candidates
cheaply (Part 1); a smallest-prime-factor sieve factors every candidate (Part 2).
Pure standard library.  Runtime: ~2-4 min.
"""
import time
from array import array
from math import isqrt
from shadow_core import liar_set, PHI24

B1 = 10**7
B2 = 2 * 10**6


def prime_sieve(B):
    sv = bytearray([1]) * (B + 1)
    sv[0] = sv[1] = 0
    for i in range(2, isqrt(B) + 1):
        if sv[i]:
            sv[i * i::i] = bytearray(len(sv[i * i::i]))
    return sv


def spf_sieve(B):
    """smallest prime factor for every n <= B."""
    spf = array('i', range(B + 1))
    for i in range(2, isqrt(B) + 1):
        if spf[i] == i:
            for j in range(i * i, B + 1, i):
                if spf[j] == j:
                    spf[j] = i
    return spf


def factor_spf(n, spf):
    f = []
    while n > 1:
        p = spf[n]
        f.append(p)
        n //= p
    return f


def legendre(a, p):
    a %= p
    return 0 if a == 0 else (1 if pow(a, (p - 1) // 2, p) == 1 else -1)


RAYS7 = (5, 7, 11, 13, 17, 19, 23)

# -------------------------------------------------- PART 1
print(f"=== PART 1: Carmichael numbers coprime to 6 up to {B1:,} ===", flush=True)
t0 = time.time()
sv = prime_sieve(B1)
carms = []
n = 25
while n <= B1:
    if n % 6 in (1, 5) and not sv[n]:           # odd, coprime to 3, composite
        if pow(2, n - 1, n) == 1:               # Fermat base-2 pseudoprime
            # factor (it is a base-2 pseudoprime: factor by trial division)
            m, f = n, []
            d = 5
            while d * d <= m:
                while m % d == 0:
                    f.append(d); m //= d
                d += 2 if d % 3 != 1 else 4      # skip multiples of 2,3 loosely
                if d % 3 == 0:
                    d += 2
            if m > 1:
                f.append(m)
            if len(f) == len(set(f)) and len(f) >= 3 and all((n - 1) % (p - 1) == 0 for p in f):
                carms.append((n, sorted(f)))
    n += 2 if n % 2 else 1

liar_counts = [len(liar_set(n)) for n, _ in carms]
witness_avg = (8 * len(carms) - sum(liar_counts)) / len(carms)
at_least_one_witness = sum(1 for c in liar_counts if c < 8)
# best Legendre agreement: per Carmichael, # of rays r with (r/p) constant over factors
best_agree, best_n = 0, None
for (n, f) in carms:
    agree = sum(1 for r in RAYS7 if len({legendre(r, p) for p in f}) == 1)
    if agree > best_agree:
        best_agree, best_n = agree, n
print(f"  count                          : {len(carms)}   (paper: 102)")
print(f"  with at least one witness      : {at_least_one_witness}/{len(carms)}  (paper: 100%)")
print(f"  average witnesses (of 8)       : {witness_avg:.2f}   (paper: 7.75)")
print(f"  best Legendre agreement        : {best_agree}/7 at n={best_n}   (paper: 4/7 at 252601)")
print(f"  elapsed                        : {time.time()-t0:.0f}s", flush=True)

# -------------------------------------------------- PART 2
print(f"\n=== PART 2: squarefree, >=3 prime factors, coprime to 6, up to {B2:,} ===", flush=True)
t0 = time.time()
spf = spf_sieve(B2)
total = with_liar = eight = 0
for n in range(5, B2 + 1):
    if n % 6 not in (1, 5):
        continue
    f = factor_spf(n, spf)
    if len(f) != len(set(f)) or len(set(f)) < 3:   # squarefree with >=3 distinct primes
        continue
    total += 1
    c = len(liar_set(n))
    if c >= 1:
        with_liar += 1
    if c == 8:
        eight += 1
print(f"  count of such composites       : {total}   (paper: 184,510)")
print(f"  with >= 1 strong liar          : {with_liar}   (paper: 130)")
print(f"  with all 8 strong liars        : {eight}   (paper: 0)")
print(f"  elapsed                        : {time.time()-t0:.0f}s", flush=True)
