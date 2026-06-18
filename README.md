# Reproducibility and verification code

Verification scripts and certificates supporting the manuscript

> **Primality testing and deterministic prime construction via modulo-6 maps and
> modulo-24 structure** — F. Bustillos, V. Leiva, C. Martin-Barreiro.

Every script here reproduces a **bounded / computational** claim of the paper —
the statements of the form *"verified up to X"*, *"exactly N up to a bound"*, or a
specific large example a referee may want to recheck. Claims that the paper proves
unconditionally are **not** reproduced here (they need no computation).

Pure Python 3 (standard library only). Run any script with `python <name>.py`.
`shadow_core.py` holds the shared primitives (shadow bases, strong-liar / Euler-liar
tests, Korselt / Carmichael test, deficiency index, factorization); the others
import it.

## Reproducibility status

**Two of the claims cannot be reproduced from this repository alone**, because
they range over the complete census of Carmichael numbers up to `1e22`, an object
that is *not regenerable* on a workstation and is taken from the literature
(R. G. E. Pinch, extended by C. Goutier; OEIS A002997, ~1.1 GB, not bundled):

- *"exactly seven eight-base-liar Carmichael numbers up to `1e22`"* (and the count
  `49,656,492`, the Pinch count `20,138,200`, the max-liars-by-range table) —
  established only by running `census_scan.py` over that external file;
- *"every Carmichael up to `1e10` is caught by GCD extraction"* — `gcd_extraction.py`
  over the same file.

`verify_8liar_carmichael.py` **does not** establish completeness: it takes the seven
numbers as given (they come out of `census_scan.py`) and verifies each is a genuine
eight-base-liar Carmichael caught by the L-check. It proves *"these seven are
eight-liars"*, not *"these are the only seven"*. The completeness rests on the cited
census. Every **other** script below is fully reproducible on a workstation and its
output matches the paper exactly.

## The eight shadow bases
For `n` coprime to 6, `sigma_r(n) = (-r * 24^{-1}) mod n` for `r` in
`{1,5,7,11,13,17,19,23}`. `n` passes the shadow MR test iff every `sigma_r` is a
strong (Miller–Rabin) liar.

## Scripts and the claim each reproduces

| Script | Reproduces (bounded claim) | Expected result | Runtime |
|---|---|---|---|
| `exhaustive_5e8.py` | Shadow test agrees with true primality for every `n` coprime to 6 up to `5e8` (no false positive). | `disagreements = 0` | ~1 h |
| `census_scan.py` † | Carmichael census to `1e22`: `49,656,492` coprime to 6; max liars `= 8` by exactly `7`; smallest `256673641562639731`; `92` seven-liars; Pinch count `20,138,200` to `1e21`; max-liars-by-range table (`1e9→4, 1e11→6, 1e13→7, 1e22→8`). | **needs external census** (not reproducible without it) | ~2 h |
| `gcd_extraction.py` † | Every Carmichael up to `1e10` is caught by GCD extraction alone (smallest 8-liar is `~2.6e17`). The extraction routine itself is checked on known Carmichaels by the suite. | **needs external census** for the full `≤1e10` sweep | seconds |
| `deficiency_index.py` | Deficiency census to `2e5`: `58` primes with `I_p=2`, `1` with `I_p=3` (`96769`); smallest `I_p=4` at `7789181`, `I_p=5` at `19249921`; `L_p=p-1` for `99.67%`. | matches paper | ~1 min |
| `group_generation.py` | Shadow bases generate `(Z/pZ)*` for every prime below the first deficient prime `6581` (**849** such primes — see *Corrected figures*); bases cover all of `Z/5Z`, `Z/7Z`; all `128` Legendre patterns realized among the `9583` primes `≤ 1e5`. | all checks pass | seconds |
| `carmichael_stats.py` | `102` Carmichael coprime to 6 up to `1e7` (avg `7.75/8` witnesses; best 3-factor Legendre agreement `4/7`, smallest at **399001** — see *Corrected figures*); `184,510` squarefree `≥3`-factor composites up to `2e6`, `130` with `≥1` liar, `0` eight-liars. | all MATCH | ~2 min |
| `korselt_triplets.py` | Korselt obstruction over the `569` primes `<2e4` with `v2(p-1)=2`: of `C(569,3)` triples, `98,806` have `≥5` depth-matches (`0` satisfy Korselt), `646,886` have exactly `4` (`5` satisfy Korselt); `340` all-QR primes `≡13 (24)` up to `5e6`, none of `C(340,3)` satisfies Korselt. | all MATCH | ~10–20 min |
| `chernick_scan.py` | Chernick family `(6k+1)(12k+1)(18k+1)`, even `k<5e7` (`n` up to `~1.6e26`): `62,734` Carmichael, max liars `6`, `0` eight-liars. | all MATCH | ~2 min |
| `semiprime_hunt.py` | Residual `8`-liar semiprimes `p=a(q-1)+1`: only the `a=2` anchor `1382809953636541 = q(2q-1)`, `q=26294581` (caught by the Δ-check); no `a≥4`. | only the anchor | ~3 min |
| `verify_8liar_carmichael.py` | The seven `8`-base-liar Carmichaels up to `1e22` (smallest `256673641562639731`) plus `psi_9`: factors, Korselt, `8/8` liars, least factor `≤ n^{1/3}`. | `ALL CHECKS PASS` | seconds |
| `verify_examples.py` | Specific examples: Euler-survivors `399001`, `488881` (`8/8` Euler liars); liar-records `4`-liar `579606301`, `6`-liar `34558584607`, `7`-liar `5588215337251`; targeted Carmichael `11239359601` (all `8` witnesses). | `ALL EXAMPLE CHECKS PASS` | seconds |
| `proth_verify.py` | The `10,006`-digit Proth prime `526623 * 2^33220 + 1`. | `N is PRIME` (witness `a=7`) | ~3 min |
| `antipodal_sieve_check.py` | Modulo-24 sieve vs Eratosthenes to `1e7`: per-ray prime counts `~83,000` (Dirichlet), exact agreement. | near-uniform, exact match | seconds |
| `section5_checks.py` | Section-5 constructive checks: Pépin–shadow, oasis-lattice twins, `kappa=1`, admissible offsets `88/345`, complement maps, inversion dichotomy. | all checks pass | seconds |

† Needs the external census data file (see below).

## Corrected figures (computational evidence)
Two figures in the manuscript were corrected to agree with the code here:
- The primes generated below the first deficient prime `6581` number **849**, not
  `1227` (`group_generation.py`). The first-failure value `6581` is correct.
- The maximum 3-factor Legendre agreement `4/7` among Carmichael `≤1e7` is attained
  at **`399001 = 31*61*211`**, not at `252601 = 41*61*101` (which aligns only `3/7`)
  (`carmichael_stats.py`).

## Census data (not bundled, ~1.1 GB)
`census_scan.py` and `gcd_extraction.py` read `carm10e22.txt` (one Carmichael number
per line, up to `1e22`). Source: R. G. E. Pinch, *The Carmichael numbers up to 1e21*
(ANTS VIII, 2007), extended to `1e22` by C. Goutier; obtainable via OEIS A002997
(`https://oeis.org/A002997`, file `a002997.7z`). Edit the `DATA` path in each script.
Integrity: the file reproduces Pinch's published count of **20,138,200** Carmichael
numbers up to `1e21` (asserted before any liar statistic is reported).

## Certificates
See `certificates/` for recorded outputs (Proth witness, the seven 8-liars, the
deficiency census, the Section-5 checks).
