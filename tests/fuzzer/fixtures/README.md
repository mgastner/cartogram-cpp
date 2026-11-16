# Saved fuzzer inputs

The fuzzer tests (`ctest -L fuzzer`, disabled in CI) run `fuzzer.py`, which writes random target CSVs for the sample maps. The four CSVs in this folder made the fuzzer tests fail on 15 Nov 2026. They are kept byte for byte, with the line endings that `fuzzer.py` writes, for follow-up investigation (PR #319).

Each CSV belongs to the sample map of the same folder name. For example:

```sh
cartogram sample_data/england_by_region_since_2023/england_by_region_since_2023.geojson \
  tests/fuzzer/fixtures/england_by_region_since_2023/england_by_region_since_2023_run3.csv
```

Successful runs with default settings:

| Input | CGAL 6.0.1, before #318 | CGAL 6.2.1, before #318 | CGAL 6.0.1, after #318 | CGAL 6.2.1, after #318 |
|---|---|---|---|---|
| England, run 3 | 1 of 11 | 0 of 11 | 8 of 8 | 3 of 3 |
| England, run 4 | 1 of 4 | 0 of 4 | 8 of 8 | 3 of 3 |
| England, run 5 | 0 of 4 | 0 of 4 | 7 of 8 | 3 of 3 |
| Malaysia, run 1 | 0 of 4 | 0 of 4 | 0 of 8 | 0 of 3 |

"Before #318" means `main` at ad0848a5 and "after #318" `main` at 1be13682, each with CGAL 6.2.1 from #319 or without it. With CGAL 6.0.1, simplification gives different outputs from run to run, so the outcome of a run varies. With CGAL 6.2.1, every run of an input gives the same output.

- **England.** Failed runs stop with "Outer boundary is not simple … called from project_with_delaunay_t". After #318, only one run of run 5 with CGAL 6.0.1 failed.
- **Malaysia.** Every run after #318 stops with "Could not converge!" after 100 integrations. Before #318, runs either stopped in the same way or with a boundary that is not simple. The values of this input range from 0.06 to 870,951, a factor of 1.5 × 10⁷.
