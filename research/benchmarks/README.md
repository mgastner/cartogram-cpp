# Benchmarks of the conformality core

## `bench_direct_evaluation.cpp` (plan step D2)

Times direct spectral evaluation at cartogram-like sizes:

- 262,144 nodes on a 256 × 256 lattice, 4 nodes per cell;
- 16², 32² and 64² modes.

For each mode set it measures:

- the modal derivatives F;
- one sweep of v and D_y v over all nodes;
- one implicit-midpoint substep of all nodes;
- one outer optimizer iteration.

The program is not part of the CMake build. To build and run it, first build
`conformality_core` in `build/Release`, then from the repository root:

```sh
c++ -std=c++20 -O3 -march=native -DNDEBUG -I include -I tests/unit \
  research/benchmarks/bench_direct_evaluation.cpp \
  build/Release/libconformality_core.a -o bench_direct_evaluation
./bench_direct_evaluation
```

### Results, 2026-11-06

Run at commit `0bd6f25a` on an Apple M2 (4 performance + 4 efficiency
cores), macOS 27.2, Apple clang 21.0.0. The program is single-threaded.
`conformality_core` was built with the repository's Release flags (`-O3
-march=native -ffp-contract=off`).

| Modes | F | v and D_y v sweep | One substep |
|---|---|---|---|
| 16² | 0.043 s | 0.060 s | 0.41 s |
| 32² | 0.137 s | 0.207 s | 1.32 s |
| 64² | 0.362 s | 0.821 s | 5.32 s |

- Throughput was 6–12 GFLOP/s, and the cost grows linearly with the number
  of modes.
- A substep needed at most 5 fixed-point iterations, so it costs about 6.5
  sweeps: the predictor, the iterations, and a final evaluation for the
  tangent.
- The manufactured map is barely distorted at this scale (E0 ≈ 1.7e-6), so
  each outer iteration took a single substep. These timings therefore do not
  predict the run time on a real cartogram. That depends on substep and
  iteration counts that only a real T0 gives.
