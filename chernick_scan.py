#!/usr/bin/env python3
"""
chernick_scan.py

Chernick Carmichael numbers  n = (6k+1)(12k+1)(18k+1)  with EVEN k (so all three
factors are = 1 mod 4, the most liar-favorable configuration), scanned for
k < 5*10^7.  This reaches n ~ 1.6*10^26.

Paper claim: there are 62,734 such Carmichael numbers; none is an eight-base
liar; the maximum shadow-liar count is 6 of 8.

A number is a Chernick Carmichael iff all three of 6k+1, 12k+1, 18k+1 are prime.
We test 6k+1 first and skip as soon as a factor is composite, so the vast
majority of k cost a single primality test.

Pure standard library.  FULL run (k<5*10^7) is long (tens of minutes); set
KBOUND smaller for a quick structural check.  Pass a bound as argv[1] to override.
"""
import sys, time
from shadow_core import is_prime, liar_set

KBOUND = int(sys.argv[1]) if len(sys.argv) > 1 else 5 * 10**7

def run(KBOUND):
    t0 = time.time()
    count = 0
    maxliar = 0
    eight = 0
    dist = {}
    k = 2
    while k < KBOUND:
        p1 = 6 * k + 1
        if is_prime(p1):
            p2 = 12 * k + 1
            if is_prime(p2):
                p3 = 18 * k + 1
                if is_prime(p3):
                    n = p1 * p2 * p3
                    count += 1
                    c = len(liar_set(n))
                    dist[c] = dist.get(c, 0) + 1
                    if c > maxliar:
                        maxliar = c
                    if c == 8:
                        eight += 1
                        print(f"  *** 8-LIAR Chernick at k={k}: n={n}", flush=True)
        k += 2
        if k % (10**7) < 2:
            print(f"  ...k~{k:,} count={count} maxliar={maxliar} t={time.time()-t0:.0f}s", flush=True)
    print(f"\nDONE k<{KBOUND:,} in {time.time()-t0:.0f}s")
    print(f"  Chernick Carmichael count : {count}   (paper at 5*10^7: 62,734)")
    print(f"  maximum shadow-liar count : {maxliar}   (paper: 6)")
    print(f"  number of 8-liars         : {eight}   (paper: 0)")
    print(f"  liar-count distribution   : {dict(sorted(dist.items()))}")

run(KBOUND)
