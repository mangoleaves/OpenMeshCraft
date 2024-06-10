# README

## Documents

> Last updated: 2024.04.21-17:50.

* [ChangeLog](./doc/Version/ChangeLog.md)

## Verified environment

> Last updated: 2024.06.10-23:15.

|        | Host         | Compiler                         |
|--------| ------------ | -------------------------------- |
|&check; | Windows10    | VS 2019 (vc142), VS 2022 (vc143) |
|&check; | Windows11    | VS 2022 (vc143)                  |
|&check; | Ubuntu 24.04 | GCC 14.0.1 X86_64                |

## external libraries

> Last updated: 2024.04.21-17:58.

* Boost
* CGAL
* Eigen3
* GMP
* GoogleTest
* MPFR
* oneTBB
* parallel-hashmap
* shewchuk-predicates

### Unix

> Last updated: 2024.06.10-23:16.

We provide gmp and mpfr libs in the "external" directory.
You can also install them on your computer and link to new ones.

* gmp
  * Install `libgmp-dev`
  * Or use gmp in external/gmp/unix
* mpfr
  * Install `libmpfr-dev`
  * Or use mpfr in external/mpfr/unix

We do not provide boost, so you need to install it on your computer and link to it.

* boost
  * Download [boost](https://www.boost.org/users/download/), required version >= 1.78.0
  * Unpack boost and run `./booststrap.sh --prefix=/usr/`
  * Run `./b2` to compile
  * Run `sudo ./b2 install` to install
