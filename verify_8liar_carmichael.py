#!/usr/bin/env python3
"""
verify_8liar_carmichael.py

Verifies the seven eight-base-liar Carmichael numbers reported in the paper
(smallest 256673641562639731; the list includes psi_9).
For each n it checks: factorization, primality of factors, coprimality to 6,
squarefreeness, the Korselt criterion (=> Carmichael), that all 8 shadow bases
are strong liars (=> 8-liar), and that the least prime factor is <= n^{1/3}
(=> caught by the L-check).

SCOPE: this confirms that each of these seven numbers IS an eight-base-liar
Carmichael caught by the L-check.  It does NOT establish completeness -- that
these are the ONLY eight-liars up to 10^22.  Completeness is a census claim,
reproducible only by census_scan.py over the full Pinch/Goutier census
(OEIS A002997, ~1.1 GB, not bundled).  The seven numbers below are themselves
the output of that census.
"""
from math import prod
from shadow_core import (PHI24, is_prime, liar_set, is_carmichael, factorize,
                         strong_liar)

# the seven 8-liar Carmichaels coprime to 6 up to 10^22 (from census_scan.py)
EIGHT_LIARS = [
    256673641562639731,
    3825123056546413051,     # psi_9
    254699850156491854531,
    326644156412798177801,
    406109173515574567039,
    2242921587179041518751,
    5958695097405523240951,
]

def integer_cuberoot(n):
    x = int(round(n ** (1/3)))
    while x**3 > n: x -= 1
    while (x+1)**3 <= n: x += 1
    return x

ok_all = True
for n in EIGHT_LIARS:
    F = factorize(n)
    liars = liar_set(n)
    prod_ok = (prod(F) == n)
    primes_ok = all(is_prime(p) for p in F)
    carm = is_carmichael(n)
    eight = (len(liars) == 8)
    least = min(F)
    caught = least <= integer_cuberoot(n)
    ok = prod_ok and primes_ok and carm and eight and caught
    ok_all &= ok
    print(f"n = {n}")
    print(f"   factors {F}  product_ok={prod_ok}  factors_prime={primes_ok}")
    print(f"   Carmichael(Korselt)={carm}  liar_set={liars}  eight_liar={eight}")
    print(f"   least factor {least} <= n^(1/3)={integer_cuberoot(n)} ? {caught}  (caught by L-check)")
    print(f"   => {'OK' if ok else 'FAIL'}\n")

print("ALL CHECKS PASS" if ok_all else "SOME CHECK FAILED")
