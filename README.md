# cartogram-cpp: Cartogram generator in C++ [![DOI](https://zenodo.org/badge/281575635.svg)](https://zenodo.org/badge/latestdoi/281575635) [![license: GPLv3](https://img.shields.io/badge/license-GPLv3-blue)](./LICENSE)

<p align="center">
    <a href="https://go-cart.io"><img src ="img/gocart_logo.svg" width="65%"></a>
</p>

This program uses the fast flow-based method developed by Michael T. Gastner, Vivien Seguy, and Pratyush More. For more information, you may refer to the following [paper](https://www.pnas.org/content/115/10/E2156):

Gastner MT, Seguy V, More P. _Fast flow-based algorithm for creating density-equalizing map projections_. Proc Natl Acad Sci USA 115(10):E2156–E2164 (2018). <https://doi.org/10.1073/pnas.0400280101>

Data produced by code in this repository are subject to the MIT license found [here](./LICENSE) and should cite the aforementioned paper by Gastner et al. (2018).

## Usage

Run the following command (replace `your-geojson-file.geojson` file with your geographic data and `your-csv-file.csv` with your visual variables file, containing target areas for each geographic region):

``` shell
cartogram your-geojson-file.geojson your-csv-file.csv
```

-   The first argument's input is a GeoJSON or JSON file, in the standard GeoJSON format.
-   The second argument's input is a `.csv` file with data about target areas.

_Note: use the `-h` flag to display more options._

The CSV file should be in the following format:

| NAME_1     | Data (e.g., Population) | Color   |
| :--------- | :---------------------- | :------ |
| Bruxelles  | 1208542                 | #e74c3c |
| Vlaanderen | 6589069                 | #f1c40f |
| Wallonie   | 3633795                 | #34495e |

-   `NAME_1` should be the same as the identifying property's name in the GeoJSON. The rows should also have the same data as is present in the identifying property.
-   `Data` contains the data you would like your cartogram to based on.
-   `Color` is the color you would like the geographic region to be. Colors may be represented in the following manner:

    1.  `cornflowerblue`: html color codes supported by `CSS3` (case-insensitive), full list of supported colors may be found in the "Extended colors" section of [web colors](https://en.wikipedia.org/wiki/Web_colors).
    2.  `"rgb(255, 0, 120)"` or `rgb(255 0 120)` or `"255, 0, 120"` or `255 0 120`: red, green and blue values out of 255.
    3.  `#e74c3c`: hex code of color, must start with `#`.

You may find sample GeoJSON (containing geographic data) and CSV (containing information about target areas, colors and other visual variables) files in the `cartogram-cpp/sample_data` directory.

To test whether whether the program was installed successfully and is working fine, you may run the following command from the repository root:

``` shell
cartogram sample_data/world_by_country_since_2022/world_by_country_since_2022.geojson sample_data/world_by_country_since_2022/world_population_by_country_2010.csv --plot_polygons --world
```

You may inspect the resultant SVG to check if everything looks as expected.

## Development

We manage dependencies with a Python virtual environment and Conan 2. The project supports GCC, Clang, and Apple Clang, all with C++20 support. Please ensure that Python 3.10 or later and a C++20-supported compiler are installed before proceeding.

Only `Debug` build commands are shown below, but the same commands can be run with `Release` build by replacing `Debug` with `Release`.

1. Create a virtual environment with the required dependencies

``` shell
virtualenv .venv && .venv/bin/pip install -U -r requirements.txt
```
and activate it:

``` shell
source .venv/bin/activate
```

2. Setup Conan

``` shell
.venv/bin/conan remote update conancenter --url=https://center2.conan.io
```

The following command will detect your system's profile and set it up for you. If you already have a profile set up, this may yeild an error, in which case you may skip this step.

``` shell
.venv/bin/conan profile detect
```

3. Install dependencies via Conan


<!-- Alternatively, we can run `export CMAKE_MINIMUM_POLICY_VERSION=3.5` before running the `conan` command to still have everything working and remove the python dependency -->

``` shell
.venv/bin/conan install . --build=missing -s build_type=Debug -s compiler.cppstd=20
```

4. Compile the project via CMake

Configure,
``` shell
.venv/bin/cmake -B build/Debug -S . -DCMAKE_TOOLCHAIN_FILE=build/Debug/generators/conan_toolchain.cmake -DCMAKE_BUILD_TYPE=Debug -DCMAKE_EXPORT_COMPILE_COMMANDS=ON
```

Build,
``` shell
.venv/bin/cmake --build build/Debug -j4
```

And, optionally, install the program globally:
``` shell
sudo .venv/bin/cmake --install build/Debug
```

### Tests

To run all the tests, execute the following command from the root directory of the repository:

``` shell
.venv/bin/ctest --test-dir build/Debug --output-on-failure
```

#### Unit Tests

To run only the unit tests:

``` shell
.venv/bin/ctest --test-dir build/Debug --output-on-failure -L unit
```

To run a specific unit test, specify the test's name. For example, to run the `test_string_to_decimal_converter.cpp` unit test, use:

``` shell
.venv/bin/ctest --test-dir build/Debug --output-on-failure test_string_to_decimal_converter.cpp
```

#### Stress Tests
This test will run all the maps in the `cartogram-cpp/sample_data` folder.

To run only the stress tests:

``` shell
.venv/bin/ctest --test-dir build/Debug --output-on-failure -L stress
```

#### Fuzzer Tests
Fuzzer tests run maps in the `cartogram-cpp/sample_data` folder with random data.

To run only the fuzzer tests:

``` shell
.venv/bin/ctest --test-dir build/Debug -L fuzzer --verbose
```
This test will take a while to finish.

Add `--verbose` to the command to see more details about the test results.

### Windows (Using WSL)

For Windows users, we recommend using our program through Windows Subsystem for Linux (WSL).

### Troubleshooting

- If you are unable to compile on the latest version of Ubuntu/macOS, please open an issue.
- If compilation suddenly stopped working for you, you may remove the `build` directory with `rm -rf build` and run the installation commands again.
- If running one of the commands starting with `.venv/bin/cmake` gives you an error, it is likely that a dependency was not installed correctly. Rerun the appropriate commands above to install the required dependencies and try again. If it still fails, make sure you have the virtual environment activated by running `source .venv/bin/activate` in your terminal, and then try again.
- If you get an error which mentions permission issues, try running the command that gave you the error with `sudo` prefixed. Alternatively, you may follow the next instruction.
- If you still get permission issues or VScode's `CMake: Install` does not work, make sure you own the relevant directories (i.e. `/usr/local/bin` and the working directory). You may assign ownership to your account with `sudo chown -R $(whoami) .`, replacing `.` with the directory of choice.

### Benchmarking

To benchmark the program, first install [hyperfine](https://github.com/sharkdp/hyperfine). You can install it using Homebrew on macOS:

```shell script
brew install hyperfine
```

Or using apt on Debian-based distributions:

```shell script
apt install hyperfine
```

Then, go to the `cartogram-cpp/tests` directory and run the following command:

```shell script
bash stress_test.sh
```

### Uninstallation

Go to the `cartogram-cpp` directory in your preferred terminal and execute the following command:

```shell script
sudo make uninstall -C build
```

Upon successful uninstallation, the following will be outputted:

    > Built target uninstall

Further, running `cartogram` should no longer work.

### Pushing changes to [go-cart.io](https://go-cart.io)

To push changes to production, please follow the the instructions on [go-cart-io/carotgram-docker](https://github.com/go-cart-io/cartogram-docker).


### Contributing

Contributions are highly encouraged! Please feel free to take a stab at any at any of the open issues and send in a pull request. If you need help getting setup or more guidance contributing, please @ any of the main contributors (@adisidev, @nihalzp, @mgastner) under any of the open issues (or after creating your own issue), and we'll be happy to guide you!

If you'd like to contribute to the project, please run our tests after you make any changes.

Maintainers, please make sure all the CI build and test checks pass and the performance comparison CI check results are expected before approving the pull request.
