#!/usr/bin/env python3
"""
Exhaustive verification of the shadow primality test to 5x10^8.

Claim (paper, Section 4.1(i)): no composite coprime to six up to 5x10^8 evades
all eight shadow bases without being caught by the Delta-check or the L-check;
equivalently, the complete shadow test agrees with true primality on every such n.

This script sieves primes to B = 5x10^8 and, for every odd n coprime to six in
[25, B], compares the shadow test verdict against the sieve. It reports the number
of disagreements (must be 0) and the running maximum shadow-liar count seen.
"""
import time
from math import gcd, isqrt

B = 5 * 10**8
PHI24 = (1, 5, 7, 11, 13, 17, 19, 23)
M = 223092870  # product of primes <= 23

def shadow_composite_caught(n, sieve_is_prime):
    """Return True iff the complete shadow test correctly classifies n.
    For a prime: all 8 bases pass and the L-check finds no factor -> PRIME.
    For a composite: caught by gcd-prefilter / perfect power / Delta / a witness / L-check."""
    # gcd pre-filter (small prime factor <= 23)
    g = gcd(n, M)
    caught_small = 1 < g < n
    # shadow MR: count witnesses (early exit at first witness)
    inv24 = pow(24, -1, n)
    d = n - 1; s = 0
    while d % 2 == 0: d //= 2; s += 1
    has_witness = False
    liars = 0
    for r in PHI24:
        a = (-r * inv24) % n
        x = pow(a, d, n)
        isl = (x == 1 or x == n - 1)
        if not isl:
            for _ in range(s - 1):
                x = x * x % n
                if x == n - 1: isl = True; break
        if isl: liars += 1
        else: has_witness = True
    # Delta-checks
    import math
    def is_square(t):
        r = math.isqrt(t); return r * r == t
    delta = is_square(8 * n + 1) or is_square(12 * n + 4)
    # L-check: trial division to n^{1/3} by d coprime to M
    nroot3 = round(n ** (1/3)) + 2
    lfound = False
    dd = 25
    while dd <= nroot3:
        if gcd(dd, M) == 1 and n % dd == 0:
            lfound = True; break
        dd += 2
    test_says_prime = (not caught_small) and (not has_witness) and (not delta) and (not lfound)
    return test_says_prime == sieve_is_prime, liars

def run():
    t0 = time.time()
    print(f"sieving to {B:,} ...", flush=True)
    sieve = bytearray([1]) * (B + 1)
    sieve[0] = sieve[1] = 0
    for i in range(2, isqrt(B) + 1):
        if sieve[i]:
            sieve[i*i::i] = bytearray(len(sieve[i*i::i]))
    print(f"sieve done ({time.time()-t0:.0f}s). scanning...", flush=True)
    mism = 0; maxliar = 0; tested = 0
    n = 25
    while n <= B:
        if n % 6 in (1, 5):  # coprime to 6
            tested += 1
            ok, liars = shadow_composite_caught(n, sieve[n] == 1)
            if not ok:
                mism += 1
                print(f"  MISMATCH at n={n} (sieve_prime={sieve[n]==1})", flush=True)
            if liars > maxliar:
                maxliar = liars
            if n % 50000000 < 2:
                print(f"  ...n~{n:,} tested={tested:,} mismatches={mism} maxliar={maxliar} t={time.time()-t0:.0f}s", flush=True)
        n += 2 if n % 2 == 1 else 1
    print(f"\nDONE to {B:,} in {time.time()-t0:.0f}s", flush=True)
    print(f"tested (coprime to 6) = {tested:,}", flush=True)
    print(f"disagreements with true primality = {mism}  (must be 0)", flush=True)
    print(f"maximum shadow-liar count <= 5x10^8 = {maxliar}", flush=True)

run()
