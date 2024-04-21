# ChangeLog

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
