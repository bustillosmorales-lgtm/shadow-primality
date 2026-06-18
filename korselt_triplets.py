#!/usr/bin/env python3
"""
korselt_triplets.py

The Korselt-obstruction census (Proposition "Korselt obstruction") and the
all-QR triple census.

PART 1.  Let P be the primes p < 20000 with nu_2(p-1) = 2.  Paper: |P| = 569.
For each unordered triple, the number of shadow rays that WOULD be strong liars
if the three primes formed a Carmichael number equals the number of rays r for
which the 2-adic depth  delta_p(sigma_r) = nu_2(ord_p(sigma_r))  is constant
across the three factors (Proposition "Depth criterion and universal cage").
Paper:
   triples with >= 5 depth-matches : 98,806  -- of these, satisfying Korselt: 0
   triples with exactly 4 matches  : 646,886 -- of these, satisfying Korselt: 5

PART 2.  Primes p = 13 (mod 24) up to 5*10^6 that are "all QR" (every ray
r in {5,7,11,13,17,19,23} is a quadratic residue).  Paper: 340 of them, and
none of the C(340,3) = 6,492,980 triples satisfies Korselt.

Pure standard library.  Runtime: ~5-15 min (Part 1 dominates).
"""
import time
from math import isqrt
from shadow_core import PHI24

POP = bytes(bin(i).count("1") for i in range(256))


def primes_upto(B):
    sv = bytearray([1]) * (B + 1)
    sv[0] = sv[1] = 0
    for i in range(2, isqrt(B) + 1):
        if sv[i]:
            sv[i * i::i] = bytearray(len(sv[i * i::i]))
    return [i for i in range(2, B + 1) if sv[i]]


def depth_vector(p):
    """( nu_2(ord_p(sigma_r)) )_{r in PHI24}, each in {0,1,2,...}."""
    inv24 = pow(24, -1, p)
    v = 0
    m = p - 1
    while m % 2 == 0:
        m //= 2
        v += 1
    out = []
    for r in PHI24:
        a = (-r * inv24) % p
        if a == 0:                # p | r (only p in {5,13} here): trivial witness
            out.append(-1)        # sentinel: never matches another factor's depth
            continue
        x = pow(a, m, p)          # order is a power of 2
        e = 0
        while x != 1:
            x = x * x % p
            e += 1
        out.append(e)
    return tuple(out)


def korselt(p, q, r):
    n1 = p * q * r - 1
    return n1 % (p - 1) == 0 and n1 % (q - 1) == 0 and n1 % (r - 1) == 0


def legendre(a, p):
    a %= p
    return 0 if a == 0 else (1 if pow(a, (p - 1) // 2, p) == 1 else -1)


# ============================================================= PART 1
print("=== PART 1: Korselt obstruction over primes p<20000 with nu2(p-1)=2 ===", flush=True)
t0 = time.time()
P = [p for p in primes_upto(20000) if p >= 5 and (p - 1) % 4 == 0 and (p - 1) % 8 != 0]
N = len(P)
print(f"  |P| = {N}   (paper: 569)", flush=True)

D = [depth_vector(p) for p in P]
# pairwise match masks: bit r set iff depth equal on ray r
pm = [[0] * N for _ in range(N)]
for i in range(N):
    Di = D[i]
    for j in range(i + 1, N):
        Dj = D[j]
        msk = 0
        for b in range(8):
            if Di[b] == Dj[b]:
                msk |= (1 << b)
        pm[i][j] = pm[j][i] = msk

print(f"  pair masks built ({time.time()-t0:.0f}s); scanning C({N},3) triples...", flush=True)
ge5 = exact4 = ge5_k = exact4_k = 0
for i in range(N):
    if i % 50 == 0:
        print(f"    ...i={i}/{N}  ge5={ge5} exact4={exact4}  t={time.time()-t0:.0f}s", flush=True)
    pmi = pm[i]
    Pi = P[i]
    for j in range(i + 1, N):
        mij = pmi[j]
        pmj = pm[j]
        Pj = P[j]
        for k in range(j + 1, N):
            c = POP[mij & pmj[k]]
            if c >= 4:
                if c >= 5:
                    ge5 += 1
                    if korselt(Pi, Pj, P[k]):
                        ge5_k += 1
                else:  # c == 4
                    exact4 += 1
                    if korselt(Pi, Pj, P[k]):
                        exact4_k += 1
print(f"  triples with >=5 depth-matches : {ge5}   (paper: 98,806)")
print(f"     of these satisfying Korselt : {ge5_k}   (paper: 0)")
print(f"  triples with exactly 4 matches : {exact4}   (paper: 646,886)")
print(f"     of these satisfying Korselt : {exact4_k}   (paper: 5)")
print(f"  elapsed                        : {time.time()-t0:.0f}s", flush=True)

# ============================================================= PART 2
print("\n=== PART 2: all-QR primes p=13 (mod 24) up to 5*10^6 ===", flush=True)
t0 = time.time()
RAYS7 = (5, 7, 11, 13, 17, 19, 23)
allqr = [p for p in primes_upto(5 * 10**6)
         if p % 24 == 13 and all(legendre(r, p) == 1 for r in RAYS7)]
print(f"  count of all-QR primes         : {len(allqr)}   (paper: 340)", flush=True)
kk = 0
M = len(allqr)
for i in range(M):
    for j in range(i + 1, M):
        for k in range(j + 1, M):
            if korselt(allqr[i], allqr[j], allqr[k]):
                kk += 1
print(f"  C({M},3) triples satisfying Korselt : {kk}   (paper: 0)")
print(f"  elapsed                        : {time.time()-t0:.0f}s", flush=True)
