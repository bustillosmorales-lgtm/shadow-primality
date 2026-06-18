# Reproducibility and verification code

Verification scripts and certificates supporting the manuscript
*"Primality testing and deterministic prime construction from six:
the modulo 24 structure"* (Bustillos, Leiva, Martin-Barreiro).

Pure Python 3 (only the standard library). Run any script with `python <name>.py`.
`shadow_core.py` holds the shared primitives (shadow bases, strong-liar test,
Korselt/Carmichael test, deficiency index); the others import it.

## The eight shadow bases
For `n` coprime to 6, `sigma_r(n) = (-r * 24^{-1}) mod n` for `r` in
`{1,5,7,11,13,17,19,23}`. `n` passes the shadow MR test iff every `sigma_r`
is a strong (Miller-Rabin) liar.

## Scripts

| Script | Verifies | Expected result | Runtime |
|---|---|---|---|
| `verify_8liar_carmichael.py` | The 7 eight-base-liar Carmichaels up to 10^22 (smallest `256673641562639731`, plus `psi_9`): factorization, Korselt, 8/8 liars, least factor `<= n^{1/3}` (caught by L-check). | `ALL CHECKS PASS` | seconds |
| `deficiency_index.py` | Definition 3.17 census to `2x10^5` (58 with `I_p=2`, 1 with `I_p=3`); `I_p` unbounded (`I_p=4` at 7789181, `I_p=5` at 19249921). | matches paper | ~1 min |
| `section5_checks.py` | Section 5 constructive claims: Pepin-shadow criterion, oasis-lattice density and twin certification, `kappa=1`, admissible offsets (density 88/345), the complement maps. | all checks pass | seconds |
| `proth_verify.py` | The 10,006-digit Proth prime `526623 * 2^33220 + 1`. | `N is PRIME` (Proth witness a=7) | ~3 min |
| `census_scan.py` | Full census of Carmichael numbers up to 10^22 for shadow liars; integrity cross-check against Pinch's `10^21` count. | see below | ~50 min |
| `exhaustive_5e8.py` | Shadow test agrees with true primality on every n coprime to 6 up to `5x10^8` (no false positives). | `disagreements = 0` | ~1 h |
| `semiprime_hunt.py` | Directed search for residual 8-liar semiprimes; finds only the `a=2` anchor (caught by the Delta-check), no `a>=4` example. | only the anchor | minutes |

## Census data (not bundled, ~1.1 GB)
`census_scan.py` reads `carm10e22.txt` (one Carmichael number per line, up to 10^22).
Source: R. G. E. Pinch, *The Carmichael numbers up to 10^21* (ANTS VIII, 2007),
extended to 10^22 by C. Goutier; obtainable via OEIS A002997
(`https://oeis.org/A002997`, file `a002997.7z`). Edit the `DATA` path in the script.

Integrity: the file reproduces Pinch's published count of **20,138,200** Carmichael
numbers up to 10^21 (the script asserts this before reporting liar statistics).

## Census result (independently reproduced)
- 49,656,492 Carmichael numbers coprime to 6 up to 10^22.
- Maximum shadow-liar count = **8**, attained by exactly **7** of them.
- Smallest 8-liar: `256673641562639731 = 24379 * 528191 * 19933079`.
- 92 seven-liars.
- Every 8-liar has at least three prime factors, hence a prime factor `<= n^{1/3}`,
  so the L-check detects all of them: the Carmichael branch is deterministic
  unconditionally.

## Certificates
See `certificates/` for the recorded outputs (Proth witness, the seven 8-liars,
the deficiency census).
