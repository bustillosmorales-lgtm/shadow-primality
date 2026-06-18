#!/usr/bin/env python3
"""
semiprime_hunt.py

Directed search for eight-base-liar semiprimes of the residual deficient form
    p = a(q-1) + 1,   n = p*q,   a = I_p >= 2,
which (for a >= 4) would have both factors > n^{1/3} and so evade the n^{1/3}
L-check while not being of the screened Delta-forms a=2 (q(2q-1)) or a=3 (q(3q-2)).

Over q up to the given bound it finds ONLY the a=2 anchor
    1382809953636541 = 26294581 * 52589161   (q(2q-1), caught by the Delta-check),
and NO a>=4 example -- consistent with the paper's "factoring frontier" framing
(these residuals, if they exist, lie beyond reach and detecting them within the
fixed-base framework requires a factor of n).
"""
import time
from math import gcd, isqrt
from shadow_core import is_prime, is_eight_liar, PHI24

PHI = set(PHI24)
A_SET = [2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]
QBOUND = 30_000_000     # raise for a deeper search (anchor appears at q=26,294,581)

def primes_upto(B):
    sv = bytearray([1]) * (B + 1); sv[0] = sv[1] = 0
    for i in range(2, isqrt(B) + 1):
        if sv[i]:
            sv[i*i::i] = bytearray(len(sv[i*i::i]))
    return [i for i in range(5, B) if sv[i] and i % 2 and i % 3]

def run():
    t0 = time.time()
    qs = primes_upto(QBOUND)
    found = []
    for q in qs:
        for a in A_SET:
            p = a * (q - 1) + 1
            if p <= q or not is_prime(p):
                continue
            n = p * q
            if n % 24 not in PHI:
                continue
            if is_eight_liar(n):
                found.append((a, q, p, n))
                tag = "Delta-form (caught at step 4)" if a in (2, 3) else "*** a>=4 (factoring frontier) ***"
                print(f"  8-LIAR semiprime: a={a} q={q} p={p} n={n}  {tag}", flush=True)
    print(f"\nDONE q<={QBOUND:,} in {time.time()-t0:.0f}s. 8-liar semiprimes found: {len(found)}")
    print(f"a>=4 (frontier) examples: {[f for f in found if f[0] >= 4]}")

if __name__ == "__main__":
    run()
