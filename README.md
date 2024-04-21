# README

## Documents

> Last updated: 2024.04.21-17:50.

* [ChangeLog](./doc/Version/ChangeLog.md)

## Verified environment

> Last updated: 2024.04.21-17:50.

|        | Host         | Compiler                         |
|--------| ------------ | -------------------------------- |
|&check; | Windows10    | VS 2019 (vc142), VS 2022 (vc143) |
|&check; | Windows11    | VS 2022 (vc143)                  |
|&cross; | Ubuntu 20.04 | GCC 11.4.0 X86_64                |

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

> Last updated: 2024.03.04-20:25.

* gmp
  * Install `libgmp-dev`
  * Or use gmp in external/gmp/unix
* mpfr
  * Install `libmpfr-dev`
  * Or use mpfr in external/mpfr/unix
* boost
  * Download [boost](https://www.boost.org/users/download/), required version >= 1.78.0
  * Unpack boost and run `./booststrap.sh --prefix=/usr/`
  * Run `./b2` to compile
  * Run `sudo ./b2 install` to install
