# Node offsets: replace R2 after the ICC 2027 submission

**Decision (author, 20 Nov 2026):** keep the R2 offsets for the ICC 2027
submission, and replace them with an established low-discrepancy sequence
when the project is revisited.

## Why

The conformality nodes sit at the offsets of the R2 sequence
(`make_r2_nodes` in `src/inset_state/conformality_adapter.cpp`). Node ℓ of a
cell is at (frac(½ + ℓ/θ), frac(½ + ℓ/θ²)) from the cell's lower-left corner,
with θ ≈ 1.3247 the real root of θ³ = θ + 1. The choice comes from the
research plan of 6 Nov (section C.1), committed in f65af8e6.

- **No peer-reviewed source introduces R2.** It comes from a blog post by
  Martin Roberts. The author flagged the blog's site as unsafe: do not link
  to it or fetch it.
- **Its uses in the literature are few.** Xiong et al. (ECCV 2020, depth
  completion) found R2 on a par with Halton(2, 3) and behind a fixed-N
  golden-ratio lattice. Slater (Quantum Inf. Process. 18:312, 2019) and
  Wolfe et al. (EGSR 2022) use it without testing it.
- **Peer-reviewed theory covers R2 only as a member of a family.** It is a
  Kronecker sequence (Weyl 1916). Its step is algebraic with 1, 1/θ and 1/θ²
  linearly independent over ℚ, so the discrepancy is O(N^(−1+ε)) (Schmidt,
  Acta Math. 125, 1970; Niederreiter 1972). The ICC draft cites exactly this
  (paper `main.tex`, Section "Discretization").

## Requirements for the replacement

1. **Deterministic,** bit-reproducible on every platform. No
   `std::uniform_real_distribution`, whose output depends on the
   implementation.
2. **No node seeded on an edge of the triangles of the first pass of T0.**
   DT0 jumps across those edges, so a node on one gets an arbitrary
   Jacobian. The edges join lattice points of the graded quadtree: cell
   sides, diagonals, and lines of slope ±2 or ±½ to hanging midpoints.
   Grid doublings happen only after at least four finished integrations
   (`InsetState::adjust_grid`, `n_finished_integrations_ > 3`), when the
   nodes have already moved. Only the seeded positions matter.
3. **An evaluation cloud** that shares no node with the optimization cloud.
4. **Even spread in each cell** for q² = 4 (optimization) and 16
   (evaluation) nodes, or a global placement with equal weights.

## Candidates checked (exact rational arithmetic, first 4096 points, 20 Nov 2026)

| Sequence | on a cell diagonal | on any edge line through lattice points (primitive directions, \|a\|, \|b\| ≤ 2) |
|---|---|---|
| Sobol' (2D, unscrambled) | 85 (the first three points already) | many |
| Halton (2, 3) | 0 | 0 |
| Halton (3, 5) | 0 | 0 |
| R2 (current) | 0 | 0 (no rational line at all) |

For all indices: a Halton(2, 3) node can lie on p x + q y = c (integers)
only if p is even and 3 divides q. A Halton(3, 5) node can lie on it only if
3 divides p and 5 divides q. No triangle edge of the graded quadtree has such
a direction.

**Recommendation:** Halton (1960), bases 2 and 3, the standard sequence.

- **Placement:** either N = q² · lx · ly points over the whole rectangle
  with weight lx · ly / N each, so the count per cell varies slightly; or q²
  consecutive indices per cell, which spreads the nodes less evenly than R2
  within a cell. Decide by measuring the agreement of the two clouds
  (currently 4.7e-4 to 3.5e-3 for the objective at β = 0.2).
- **Alternatives:** Owen-scrambled Sobol' with a fixed seed (more code), or
  jittered strata with a specified generator (`std::mt19937_64`, raw bits).

## Work when the change is made (estimate of 20 Nov 2026)

1. **Code,** about 1 hour: a function with the signature of `make_r2_nodes`
   (called twice, `conformality_adapter.cpp:622–623`). A unit test: no node
   on an edge line, disjoint clouds, determinism. A comment at the generator
   naming this note. Unit suite; flag-off fixtures unchanged.
2. **Reruns,** one night: everything numeric depends on the nodes.
   - paper run set at β = 0.2 with mesh checks and figures, and the
     β = 1 set;
   - timing repeats, first, on a cool, idle machine;
   - serial check;
   - the 26-map rollouts at 32 modes for β = 1, 0.5, 0.2 and 0.1;
   - the 64-mode rollout and its diagnosis, if still reported;
   - the concept figure.
3. **Regeneration,** about 3 hours: tables, map and China figures,
   convergence figure, text numbers (`text_numbers_e_half.py`,
   `region_turns.py`, `summarize_pop.py`), and the "Discretization" text.
4. **Expect changes:** the optimizer's path depends on the nodes, so every
   number moves a little. Recheck that all 26 maps still pass at 32 modes.

Illustration of the current clouds and of discrepancy: paper repository,
`notes/node_clouds.html`.
