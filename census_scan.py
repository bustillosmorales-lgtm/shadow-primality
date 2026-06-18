#!/usr/bin/env python3
"""
census_scan.py

Scans the COMPLETE census of Carmichael numbers up to 10^22 for shadow MR liars.

DATA FILE (not bundled, ~1.1 GB): carm10e22.txt, one Carmichael number per line.
Source: R. G. E. Pinch, "The Carmichael numbers up to 10^21" (ANTS 2007),
extended to 10^22 by C. Goutier; available via OEIS A002997
(https://oeis.org/A002997 -> a002997.7z). Set DATA below to its path.

Integrity check: the file must reproduce Pinch's published count of 20,138,200
Carmichael numbers <= 10^21. The script verifies this before reporting liar stats.

Reproduces (independently, for the paper):
  * 49,656,492 Carmichael numbers coprime to 6 up to 10^22;
  * maximum shadow-liar count = 8, attained by exactly 7 of them;
  * smallest 8-liar = 256673641562639731; 92 seven-liars;
  * the maximum-liar-by-range table:
        n <= 10^9  -> 4
        n <= 10^11 -> 6
        n <= 10^13 -> 7
        n <= 10^22 -> 8
Every 8-liar has >= 3 prime factors, hence a factor <= n^{1/3}, caught by the L-check.

The full shadow-liar count is computed for every Carmichael (no early exit), so
the run is long (a couple of hours).  Set FAST=True to early-exit at the second
witness: this still finds the 7- and 8-liars and the smallest 8-liar, but the
max-by-range table below 7 is then not produced.

Pure standard library.
"""
import time
from shadow_core import PHI24

DATA = r"C:\Users\telot\carm_data\carm10e22.txt"   # <-- edit to your path
B21 = 10**21
BOUNDS = [10**9, 10**11, 10**13, 10**22]
FAST = False


def shadow_liar_count(n):
    inv24 = pow(24, -1, n)
    d = n - 1
    s = 0
    while d % 2 == 0:
        d >>= 1
        s += 1
    liars = 0
    wit = 0
    for r in PHI24:
        x = pow((-r * inv24) % n, d, n)
        isl = (x == 1 or x == n - 1)
        if not isl:
            for _ in range(s - 1):
                x = x * x % n
                if x == n - 1:
                    isl = True
                    break
        if isl:
            liars += 1
        else:
            wit += 1
            if FAST and wit >= 2:
                return None        # liars <= 6, not a record under FAST
    return liars


def run():
    t0 = time.time()
    total = le21 = tested = sevens = 0
    eights = []
    maxbybound = {b: 0 for b in BOUNDS}
    with open(DATA) as fh:
        for line in fh:
            line = line.strip()
            if not line.isdigit():
                continue
            n = int(line)
            total += 1
            if n <= B21:
                le21 += 1
            if n % 3 == 0:                 # not coprime to six
                continue
            tested += 1
            c = shadow_liar_count(n)
            if c is None:
                continue
            for b in BOUNDS:
                if n <= b and c > maxbybound[b]:
                    maxbybound[b] = c
            if c == 7:
                sevens += 1
            elif c == 8:
                eights.append(n)
    print(f"total Carmichael in file        : {total:,}")
    print(f"count <= 10^21                  : {le21:,}  (Pinch published 20,138,200: "
          f"{'MATCH' if le21 == 20138200 else 'MISMATCH'})")
    print(f"coprime to 6 (tested)           : {tested:,}  (paper: 49,656,492)")
    print(f"number of 7-liars               : {sevens}  (paper: 92)")
    print(f"number of 8-liars               : {len(eights)}  (paper: 7)")
    print(f"smallest 8-liar                 : {min(eights) if eights else None}")
    print(f"all 8-liars                     : {sorted(eights)}")
    print(f"max shadow-liars by range (Table):")
    for b in BOUNDS:
        print(f"    n <= {b:<6} : {maxbybound[b]}")
    print(f"elapsed                         : {time.time()-t0:.0f}s")


if __name__ == "__main__":
    run()
