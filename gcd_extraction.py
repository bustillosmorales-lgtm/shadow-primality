#!/usr/bin/env python3
"""
gcd_extraction.py

Claim (Section 3, "Closing the Carmichael gap" / Remark on shadow-Lucas):
every Carmichael number up to 10^10 is caught by GCD extraction alone -- i.e.
at least one shadow base is a strong witness and, during its Miller-Rabin
squaring chain, some intermediate value x = sigma_r^{2^j d} yields a nontrivial
gcd(x-1, n) or gcd(x+1, n) (Theorem "GCD extraction").  This holds because the
smallest eight-base liar Carmichael number is ~2.6*10^17, far above 10^10, so
every Carmichael below 10^10 has at least one witness.

DATA FILE (not bundled, ~1.1 GB): the same carm10e22.txt used by census_scan.py
(Pinch/Goutier, OEIS A002997).  We stream it, keep the Carmichael numbers <= 10^10
that are coprime to six, and run GCD extraction on each.

Pure standard library.  Runtime: seconds once the (small) <=10^10 slice is read.
"""
from math import gcd
from shadow_core import PHI24, factorize

DATA = r"C:\Users\telot\carm_data\carm10e22.txt"   # <-- same file as census_scan.py
B = 10**10


def gcd_extract(n):
    """Try to recover a nontrivial factor of Carmichael n via the shadow bases.
    Returns a factor, or None if every base is a strong liar (no extraction)."""
    inv24 = pow(24, -1, n)
    d = n - 1
    s = 0
    while d % 2 == 0:
        d //= 2
        s += 1
    for r in PHI24:
        a = (-r * inv24) % n
        x = pow(a, d, n)
        if x in (1, n - 1):
            continue                      # this base is a liar so far
        g = gcd(x - 1, n)
        if 1 < g < n:
            return g
        for _ in range(s - 1):
            x = x * x % n
            g = gcd(x - 1, n)
            if 1 < g < n:
                return g
            g = gcd(x + 1, n)
            if 1 < g < n:
                return g
            if x == n - 1:
                break
    return None


def run():
    tested = caught = 0
    failures = []
    with open(DATA) as fh:
        for line in fh:
            line = line.strip()
            if not line.isdigit():
                continue
            n = int(line)
            if n > B:
                continue                  # file is sorted; could 'break' if guaranteed sorted
            if n % 3 == 0:
                continue                  # not coprime to six
            tested += 1
            f = gcd_extract(n)
            if f is not None and n % f == 0:
                caught += 1
            else:
                failures.append(n)
    print(f"Carmichael coprime to 6 with n <= 10^10 tested : {tested}")
    print(f"caught by GCD extraction                       : {caught}")
    print(f"not caught (should be 0 below 10^10)           : {len(failures)}  {failures[:10]}")
    print("RESULT:", "ALL caught by GCD extraction" if not failures else "SOME NOT caught")


if __name__ == "__main__":
    run()
