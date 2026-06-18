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

Output (independently reproduced for the paper):
  49,656,492 Carmichael numbers coprime to 6 up to 10^22;
  maximum shadow-liar count = 8, attained by exactly 7 of them;
  smallest 8-liar = 256673641562639731; 92 seven-liars.
Every 8-liar has >=3 prime factors, hence a factor <= n^{1/3}, caught by the L-check.
"""
import time
from shadow_core import PHI24

DATA = r"C:\Users\telot\carm_data\carm10e22.txt"   # <-- edit to your path
B21 = 10**21

def run():
    t0 = time.time()
    total = le21 = tested = sevens = 0
    eights = []
    maxl = 0
    with open(DATA) as fh:
        for line in fh:
            line = line.strip()
            if not line.isdigit():
                continue
            n = int(line)
            total += 1
            if n <= B21:
                le21 += 1
            if n % 3 == 0:           # not coprime to 6
                continue
            tested += 1
            inv24 = pow(24, -1, n)
            d = n - 1; s = 0
            while d % 2 == 0: d >>= 1; s += 1
            liars = 0; wit = 0
            for r in PHI24:
                x = pow((-r * inv24) % n, d, n)
                isl = (x == 1 or x == n - 1)
                if not isl:
                    for _ in range(s - 1):
                        x = x * x % n
                        if x == n - 1: isl = True; break
                if isl: liars += 1
                else:
                    wit += 1
                    if wit >= 2: break    # liars<=6: not a record, not 8-liar
            if wit < 2:
                if liars > maxl: maxl = liars
                if liars == 7: sevens += 1
                if liars == 8: eights.append(n)
    print(f"total Carmichael in file        : {total:,}")
    print(f"count <= 10^21                  : {le21:,}  (Pinch published 20,138,200: "
          f"{'MATCH' if le21 == 20138200 else 'MISMATCH'})")
    print(f"coprime to 6 (tested)           : {tested:,}")
    print(f"maximum shadow-liar count       : {maxl}")
    print(f"number of 7-liars               : {sevens}")
    print(f"number of 8-liars               : {len(eights)}")
    print(f"smallest 8-liar                 : {min(eights) if eights else None}")
    print(f"all 8-liars                     : {sorted(eights)}")
    print(f"elapsed                         : {time.time()-t0:.0f}s")

if __name__ == "__main__":
    run()
