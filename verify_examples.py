#!/usr/bin/env python3
"""
verify_examples.py

Verifies every SPECIFIC numerical example the paper exhibits (each is a concrete
number a referee can recheck), independent of any large scan:

  * Euler-survivors 399001 = 31*61*211 and 488881 = 37*73*181:
    Carmichael, and all eight shadow bases are Euler liars (a^((n-1)/2) = +-1),
    which is why a Fermat/Euler-type shadow criterion does not catch them and the
    strong-MR refinement is needed.  We also report the strong-MR liar count.

  * Liar-count records (Table "maximum shadow liars"):
      4-liar  579606301      = 109 * 541 * 9829
      6-liar  34558584607    = 443 * 1327 * 58787   (all factors = 3 mod 4)
      7-liar  5588215337251  = 3571 * 20231 * 77351

  * Targeted Carmichael 11239359601 = 281 * 4201 * 9521: the three factors share
    the Legendre pattern (1,1,-1,-1,1,-1,-1) over r in {5,7,11,13,17,19,23} but
    the ray (-24/p_i) differs, so condition (i) of the alignment theorem fails and
    all eight bases are witnesses (liar count 0).

Pure standard library.  Runtime: seconds.
"""
from math import prod
from shadow_core import (PHI24, is_prime, is_carmichael, factorize,
                         liar_set, shadow_bases)


def euler_liar(a, n):
    """a is an Euler liar for n iff a^((n-1)/2) = +-1 (mod n)."""
    x = pow(a, (n - 1) // 2, n)
    return x == 1 or x == n - 1


def check(n, expected_factors, expected_liars=None, label=""):
    F = sorted(factorize(n))
    fac_ok = (F == sorted(expected_factors)) and prod(F) == n
    carm = is_carmichael(n)
    L = liar_set(n)
    line = (f"  n = {n}  {label}\n"
            f"     factors {F}  matches stated: {fac_ok}  Carmichael: {carm}\n"
            f"     strong-MR liar rays: {L}  (count {len(L)})")
    ok = fac_ok and carm
    if expected_liars is not None:
        ok = ok and (len(L) == expected_liars)
        line += f"  expected {expected_liars}: {len(L) == expected_liars}"
    print(line)
    return ok

ok_all = True

print("=== Euler-survivors (all 8 shadow bases are Euler liars) ===")
for n, fac in [(399001, [31, 61, 211]), (488881, [37, 73, 181])]:
    sb = shadow_bases(n)
    euler = sum(euler_liar(a, n) for a in sb)
    L = liar_set(n)
    carm = is_carmichael(n)
    print(f"  n = {n} = {'*'.join(map(str, fac))}  Carmichael: {carm}")
    print(f"     Euler liars among the 8 shadow bases: {euler}/8  (paper: survives every shadow base)")
    print(f"     strong-MR liar rays: {L}  (count {len(L)} -> MR refinement catches it: {len(L) < 8})")
    ok_all &= carm and euler == 8

print("\n=== liar-count records ===")
ok_all &= check(579606301,    [109, 541, 9829],    4, "(4-liar)")
ok_all &= check(34558584607,  [443, 1327, 58787],  6, "(6-liar, all factors = 3 mod 4)")
ok_all &= check(5588215337251,[3571, 20231, 77351],7, "(7-liar)")

print("\n=== targeted Carmichael (alignment condition (i) fails -> all 8 witnesses) ===")
n = 11239359601
F = sorted(factorize(n))
RAYS7 = (5, 7, 11, 13, 17, 19, 23)
def leg(a, p):
    a %= p
    return 0 if a == 0 else (1 if pow(a, (p - 1) // 2, p) == 1 else -1)
print(f"  n = {n}  factors {F}  Carmichael: {is_carmichael(n)}")
for p in F:
    patt = tuple(leg(r, p) for r in RAYS7)
    ray = leg(-24, p)   # (-24/p) = +1 primary ray, -1 antipodal
    print(f"     p={p}: (r/p) for r in {RAYS7} = {patt}   (-24/p) = {ray:+d}")
L = liar_set(n)
print(f"     strong-MR liar rays: {L}  (all 8 are witnesses: {len(L) == 0})")
ok_all &= is_carmichael(n) and len(L) == 0

print("\n" + ("ALL EXAMPLE CHECKS PASS" if ok_all else "SOME EXAMPLE CHECK FAILED"))
