# Baseline runs, 2026-11-06

These are the reference outputs of the unmodified cartogram-cpp, recorded
before any conformality work. They correspond to step D0 of the research
plan.

## Code and toolchain

- **Commit:** `6d76b8bd`, branch `research/conformality`. The code is the
  same as `9200b14b`; `6d76b8bd` only adds one `.gitignore` line.
- **Binary:** `build/Release/cartogram`, version string 25.10.
  - Configured by
    `.venv/bin/cmake -B build/Release -S . -DCMAKE_TOOLCHAIN_FILE=build/Release/generators/conan_toolchain.cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_EXPORT_COMPILE_COMMANDS=ON`.
  - The existing configure step was re-run because the previous one pointed
    to the removed `MacOSX26.5.sdk`.
- **Compiler:** Apple clang 21.0.0 (clang-2100.4.1.1), Xcode 27.2 (27B52),
  SDK `MacOSX27.2.sdk`.
- **Build tools:** CMake 3.30.0, Conan 2.16.1.
  - The Conan profile says apple-clang 17 and cppstd gnu17. That only affects
    the package IDs: the actual compile uses `-std=c++20`.
- **Dependencies** (all from the Conan cache, not Homebrew):
  - CGAL 6.0.1
  - Boost 1.83 (header-only)
  - FFTW 3.3.10 (static)
  - Eigen 3.4.0 (pulled in by CGAL)
  - nlohmann_json 3.11.3
  - argparse 3.2
- **Machine:** Apple M2 (arm64), macOS 27.2.
- **Unit tests:** `ctest -L unit` passed 6 of 6.

## Runs

Every run was launched from its own subdirectory with
`build/Release/cartogram <geojson> <csv> [flags] --export_time_report`.

| Directory | Flags | Integrations | Grid | Final max area error | Time | Deterministic? |
|---|---|---|---|---|---|---|
| `conterminous_usa/` | defaults | 6 | 256² | 0.99 % (New Hampshire) | 0.29 s | **no** |
| `china/` | defaults | 9 | 256² → 512² after pass 7 | 0.33 % (Macao) | 0.51 s | **no** |
| `conterminous_usa_nosimplify/` | `-P 100000000` | 6 | 256² | 0.77 % (Oregon) | 0.27 s | yes (8 of 8) |
| `china_nosimplify/` | `-P 100000000` | 8 | 256², no doubling | 0.55 % (Inner Mongolia) | 0.52 s | yes (3 of 3) |
| `conterminous_usa_nosimplify_n128/` | `-P 100000000 -n 128` | 11 | 128² → 256² (one doubling) | 0.38 % (Oregon) | 0.40 s | yes (3 of 3) |

**Research runs use `-P 100000000`, so simplification is off** (decided
2026-11-06).

- Every pass rasterizes the density from the current polygons. As a result,
  simplification changes T0 itself, not only the drawn output.
- Simplification is also what makes T0 non-reproducible.
- With simplification off, China no longer triggers `adjust_grid`. The run
  `conterminous_usa_nosimplify_n128/` is the real-data fixture for the grid
  doubling (×2) path instead.

Data sources:

- USA: `sample_data/conterminous_usa_by_state_since_1959/`, with
  `conterminous_usa_population_by_state_2020.csv`.
- China: `sample_data/china_by_province_with_chinese_taipei_since_1998/`,
  with `china_population_by_province_2020.csv`.

## SHA-256 of the archived outputs

```
d43e9549ce8768dec86caa28247020964d7e6e1819bfa7e298311a74952df1a2  conterminous_usa/conterminous_usa_population_by_state_2020_cartogram.geojson
1bcf92e079b34410152b3866f8ea5bf01dbf503cc6d73d9b65b82caa12562e2f  china/china_population_by_province_2020_cartogram.geojson
b8391ad52d490739da80df0c173db4309e35ca8165c14148464b3a1e3924cf95  conterminous_usa_nosimplify/conterminous_usa_population_by_state_2020_cartogram.geojson
d84f8fc3e5373edf403146447a12ac16884455c1635a0c7d8fe751055b4b8ffc  china_nosimplify/china_population_by_province_2020_cartogram.geojson
091efa6ad41f7e8e8bc8aa9f5a23b90ab3c173f702a3d6242ea263c2f61fef30  conterminous_usa_nosimplify_n128/conterminous_usa_population_by_state_2020_cartogram.geojson
```

## Default settings are not reproducible run to run

**What was observed:**

- **USA with default settings:** 10 runs gave 6 distinct outputs.
  - Every run exited with code 0.
  - Coordinates differed by up to 2.3·10⁻⁵ lattice units on a 256 grid.
  - The printed area errors agree to about 6 digits.
- **China with default settings:** 5 runs gave 4 distinct outputs.
- **`--export_time_report` is not the cause.** Runs with and without it hit
  the same set of outputs.
- **Uninitialized heap reads are not the cause.** Setting
  `MallocPreScribble=1` did not make the output reproducible.

**What isolates it:**

| Setting | Runs | Distinct outputs |
|---|---|---|
| `-S` (no simplification, no densification) | 8 | 1 |
| `-P 100000000` (densification on, simplification never triggered) | 8 | 1 |

The source is therefore `InsetState::simplify`
(`src/inset_state/simplify_inset.cpp`), which uses CGAL
`Polyline_simplification_2` on a `Constrained_triangulation_plus_2`.

- **Suspected mechanism (not verified):** ties are broken by vertex-handle
  address, which changes under ASLR.

**Consequence:** with default settings, the archived `conterminous_usa/` and
`china/` outputs are one sample each from a small set of possible outputs.
They cannot serve as byte-identical regression references. The
`*_nosimplify` runs can.
