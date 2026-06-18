# shadow-primality

Reference code and library for the article

**Primality testing and deterministic prime construction via modulo-6 maps and modulo-24 structure**
F. Bustillos, V. Leiva, and C. Martin-Barreiro.

## Overview

A single arithmetic structure built from the seed `6 = 2 x 3`, read in two complementary directions:

- **Forward (testing).** The three quadratic maps of six and their modulo-24 refinement yield a
  shadow function in closed form. Its eight values act as Miller-Rabin witness bases derived from the
  candidate itself. The resulting test is deterministic for every prime, every prime power, every
  Carmichael number, and every semiprime of deficiency index at most three, leaving only a class of
  balanced semiprimes at the factoring frontier.
- **Reverse (construction).** The same maps locate primes by position. A Chinese-remainder-theorem
  oasis-lattice constructor generates the primes of any interval, including twin primes, by modular
  arithmetic alone and certifies them by construction. A Pepin-type criterion certifies the
  generalized Fermat family of base `576 = 24^2` from a single exponentiation.

## Contents

The Python code used for the computational study (Python 3.12, built-in arbitrary-precision integers)
is provided here.

## Citation

If you use this code, please cite the article.
