# cartogram-cpp: Cartogram generator in C++ [![DOI](https://zenodo.org/badge/281575635.svg)](https://zenodo.org/badge/latestdoi/281575635) [![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](./LICENSE)

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

We manage dependencies with a Python virtual environment and Conan 2. The project uses Clang with C++20 support. Please ensure that Python 3.10 or later and a Clang++ (preferably Clang 20) compiler with C++20 support are installed before proceeding.

Only `Debug` build commands are shown below, but the same commands can be run with `Release` build by replacing `Debug` with `Release`.

### Linux and macOS

#### Create a virtual environment and activate it

```
python3 -m venv .venv
source .venv/bin/activate
```

#### Install dependencies while in the virtual environment

```
pip install --upgrade pip wheel conan==2.16.1 cmake==3.30.0
```

#### Setup Conan

```
conan remote update conancenter --url=https://center2.conan.io
conan profile detect
```

#### Install dependencies via Conan


<!-- Alternatively, we can run `export CMAKE_MINIMUM_POLICY_VERSION=3.5` before running the `conan` command to still have everything working and remove the python dependency -->

```
conan install . --build=missing -s build_type=Debug -s compiler.cppstd=20
```

#### Configure the project

```
.venv/bin/cmake -B build/Debug -S . -DCMAKE_TOOLCHAIN_FILE=build/Debug/generators/conan_toolchain.cmake -DCMAKE_BUILD_TYPE=Debug -DCMAKE_EXPORT_COMPILE_COMMANDS=ON
```

#### Build the project

```
.venv/bin/cmake --build build/Debug -j4
```

#### Install the project (optional)

```
sudo .venv/bin/cmake --install build/Debug

```

### Tests

To run all the tests, execute the following command from the root directory of the repository:

```
ctest --test-dir build/Debug --output-on-failure
```

#### Unit Tests

To run only the unit tests:

```
ctest --test-dir build/Debug --output-on-failure -L unit
```

```
ctest --test-dir build/Debug --output-on-failure test_string_to_decimal_converter.cpp
```

#### Stress Tests
This test will run all the maps in the `cartogram-cpp/sample_data` folder.

To run only the stress tests:

```
ctest --test-dir build/Debug --output-on-failure -L stress
```

#### Fuzzer Tests
Fuzzer tests run maps in the `cartogram-cpp/sample_data` folder with random data.

To run only the fuzzer tests:

```
ctest --test-dir build/Debug -L fuzzer --verbose
```
This test will take a while to finish.

Add `--verbose` to the command to see more details about the test results.

### Windows (Using WSL)

For Windows users, we recommend using our program through Windows Subsystem for Linux (WSL).


#### Installing using Docker

If you prefer, you may run the program in an isolated environment using [Docker Desktop](https://www.docker.com/products/docker-desktop/).

1. Start the Docker container by executing the command below. Please note, the first time you run this command may take longer as Docker will download the necessary images and build the container. Subsequent runs should get the container up and running right away.

```bash
docker compose up -d
```

2. Once the Docker container is up and running, you may access the container's shell by executing:

```bash
docker exec -it cartogram-cpp /bin/bash
```

From there, you can run commands with `cartogram` as usual. By default, you should find yourself in the `/cartogram/output` directory, which is mounted to the `output` directory in the repository root on your local environment. This means that any files you generate in `/cartogram/output` should also be available in the `output` directory in the repository root, and vice-versa. You can test out how this works by running the following command:

```bash
cartogram ../sample_data/world_by_country_since_2022/world_by_country_since_2022.geojson ../sample_data/world_by_country_since_2022/world_population_by_country_2010.csv --plot_polygons --world
```

To stop the Docker container, execute:

```bash
docker compose down
```

To compile code changes within the Docker container, run the following command from the `/cartogram` directory:

```bash
bash build.sh
```

### Troubleshooting

- If you are unable to copmile on the latest version of Ubuntu, please open an issue. In the meanwhile, follow the instructions for installation via Docker.
- If compilation suddenly stopped working for you, you may remove the `build` directory with `rm -rf build` and run the installation commands again.
- If running `cmake -B build` gives you an error, it is likely that a dependency was not installed correctly. Rerun the appropriate commands above to install the required dependencies and try again.
- If you get an error which mentions permission issues, try running the command that gave you the error with `sudo` prefixed, as done with `sudo make install -C build` above.
- If `cmake` complains that it could not find a particular library, please try uninstalling it and installing it again. After reinstalling it, please also unlink it and link it with the `--force` flag.
- If you get errors related to CGAL, it's likely you have another version of CGAL installed on your computer that is getting chosen instead of the one contained as a submodule within this repository. It's also possible that when cloning this repository, the `--recurse-submodule` flag was missing. Try running `git submodule init` and `git submodule update` in the root directory of the repository.
- If VScode's `CMake: Install` does not work, make sure you own `/usr/local/bin` and the working directory. You may assign ownership to your account with `sudo chown -R $(whoami) .`, replacing `.` with the directory of choice.

### Testing

If you'd like to contribute to the project, please run our tests after you make any changes.

To run the unit tests, execute the following command:

```shell script
ctest --verbose
```

To learn more about the tests, you may go to the `cartogram-cpp/tests` directory and read the `README.md` file.

Additionally, you may go to the `cartogram-cpp/tests` directory and run the following command:

```shell script
bash stress_test.sh
```

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

Maintainers, please make sure to run the "Build and Release" workflow under GitHub Actions before approving the pull request. You may delete the newly created release before merging the pull-request. Another release should be automatically created after merging with main.
