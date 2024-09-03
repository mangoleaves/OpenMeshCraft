# README

## Main projects

> Last updated: 2024.09.02-15:32.

### Mesh arrangements

* Reference paper
  * Jia-Peng Guo and Xiao-Ming Fu. 2024. Exact and Efficient Intersection Resolution for Mesh Arrangements. *ACM Trans. Graph*. 43, 6 (December 2024), 14 pages.
  * Cherchi, G., Livesu, M., Scateni, R. and Attene, M. Fast and robust mesh arrangements using floating-point arithmetic. ACM Transactions on Graphics, 39, 6 (2020), 1-16.
* Source codes
  * `src/OpenMeshCraft/Arrangements`
* Running examples
  * `test/Executables/arrangements.cpp` (in CMake target `OpenMeshCraft-Arrangements`)
  * `test/Arrangements/test_arrangements.cpp` (in CMake target `OpenMeshCraft-Test`)
* Tested on Thingi10k dataset and one dataset from the first reference paper.

### Mesh boolean

* Reference paper
  * Cherchi, G., Pellacini, F., Attene, M. and Livesu, M. Interactive and Robust Mesh Booleans. *ACM Transactions on Graphics*, 41, 6 (2022), 1-14.
* Source codes
  * `src/OpenMeshCraft/Boolean`
* Running examples
  * `test/Boolean/test_boolean.cpp` (in CMake target `OpenMeshCraft-Test`)
* Haven't been tested thoroughly.

## Documents

> Last updated: 2024.04.21-17:50.

* [ChangeLog](./ChangeLog.md)

## Build

> Last updated: 2024.09.02-15:39.

1. Configure the path to the Boost library. All other third-party libraries are already included.
2. Run the CMake configuration and build the target you intend to run.

### Verified environment

> Last updated: 2024.06.10-23:15.

|        | Host         | Compiler                         |
|--------| ------------ | -------------------------------- |
|&check; | Windows10    | VS 2019 (vc142), VS 2022 (vc143) |
|&check; | Windows11    | VS 2022 (vc143)                  |
|&check; | Ubuntu 24.04 | GCC 14.0.1 X86_64                |

### external libraries

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

## License

This project is licensed under the terms of the GNU General Public License (GPL) and the GNU Lesser General Public License (LGPL).

* The source code and related materials are licensed under the [GNU General Public License v3.0](https://www.gnu.org/licenses/gpl-3.0.html). For more information, see the LICENSE file.
* Some components of this project are licensed under the [GNU Lesser General Public License v3.0](https://www.gnu.org/licenses/lgpl-3.0.html), which allows for linking with proprietary modules.

By contributing to this project, you agree that your contributions will be licensed under the same terms as the rest of the project.

## Bug report

Please report bugs by opening a GitHub issue and providing the model that caused the error. While I strive to minimize bugs, some may have gone untested during development, and I will do my best to correct them `:)`.
