# ChangeLog

## 1.0.2 (2024.09.02)

Algorithm Improvements:

- Improve implementations of multi-precision arithmetics.
- Compress expansions and update interval numbers by expansions to improve efficiency.
- Optimize implicit point order to improve accuray of semi-static filter.

Docs:

- Update README.md and test/README.md

Bug Fix:

- Fix bugs in *Mesh arrangements* and *Mesh boolean*

## 1.0.1 (2024.06.10)

New Features:

- No

Algorithm Improvements:

- Add more controls for *Mesh Arrangements* to enable/disable features.
- Add more profiling codes and tools for *Mesh Arrangements*. Parallel profiling is available.
- Improve efficiency of *Mesh Arrangements*:
  - Cache SSF results in implicit points.
  - Switch expansion to exact rational numbers if expansion is too complex.
  - Use auxiliary seprator instead of segment itself to locate the segment. It simplify orientOn2D with at most three implicit points to Orient3D with at most one implicit point.
  - Cut ears containing TPI points first, then cut remaining ears. It avoids OrientOn2D with more than two TPI points.
  - Update paramters of adaptive OcTree
  - Remove useless features

Bug Fix:

- Fix bugs in *Mesh Arrangements*.

Miscellaneous:

- Support compiling on ubuntu-24.04 with gcc14.

## 1.0.0 (2024.04.21)

New Features:

- Geometric primitives and predicates (needs to be refactored for easier use)
  - Indirect predicates
  - Offset indirect predicates
- Predicates generator
  - Automatically generate indirect predicates and offset indirect predicates
- Geometric kernels
  - Approximate predicates and approximate constructions (not fully tested)
  - Exact (offset) indirect predicates and approximate constructions (fully tested in mesh arrangements)
- Intersection tests for basic geometry primitives (incomplete for all known primitives)
- BVHs
  - AABB Tree
  - Kd Tree
  - Binary Tree
  - Orthogonal Tree
  - Adaptive Orthogonal Tree
- Number types used for adaptive precision computation
  - Interval number (from IndirectPredicates, Macro Attene)
  - Rational number (wrap of gmp, mpfr, boost::multiprecision)
  - Expansion floating point number (from IndirectPredicates, Macro Attene)
  - Big numbers (from IndirectPredicates, Macro Attene)
- Mesh arrangements (stable, tested on Thingi10k dataset and a complex dataset)
- Mesh boolean (unstable, not fully tested)
