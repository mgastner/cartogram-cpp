# cartogram-cpp: Cartogram generator in C++ [![DOI](https://zenodo.org/badge/281575635.svg)](https://zenodo.org/badge/latestdoi/281575635) [![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](./LICENSE)

<p align="center">
    <a href="https://go-cart.io"><img src ="img/gocart_logo.svg" width="65%"></a>
</p>

This program uses the fast flow-based method developed by Michael T. Gastner, Vivien Seguy, and Pratyush More. For more information, you may refer to the following [paper](https://www.pnas.org/content/115/10/E2156):

Gastner MT, Seguy V, More P. _Fast flow-based algorithm for creating density-equalizing map projections_. Proc Natl Acad Sci USA 115(10):E2156–E2164 (2018). <https://doi.org/10.1073/pnas.0400280101>

Data produced by code in this repository are subject to the MIT license found [here](./LICENSE) and should cite the aforementioned paper by Gastner et al. (2018).

While cloning this repository, please ensure you use the `--recurse-submodules` flag like so:
-
    git clone --recurse-submodules https://github.com/mgastner/cartogram-cpp.git

## Dependencies

Please note, we only support UNIX-based systems, and have only tested on macOS, Linux, and GNU.

### macOS

#### Installing Homebrew

Install [homebrew](brew.sh) by running the following command:

    /bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"

#### Installing dependencies through Homebrew

Install pkg-config, boost, fftw, nlohmann-json, and cmake by running the following command:

    brew install libomp pkg-config boost fftw nlohmann-json cmake cairo

### Debian-based distributions (Ubuntu, Arch Linux etc.)

#### Installing relevant dependencies through apt:

Have a look through to apt-requirements.txt if you'd like to see what all will be installed. Then, run the following commands to install all dependencies through apt:

    apt install -y g++-11 build-essential cmake libboost-all-dev nlohmann-json3-dev libomp-dev libfftw3-dev libcairo2-dev

### Installation

Go to the `cartogram-cpp` directory in your preferred terminal and execute the following commands.

    cmake -B build
    make -C build
    sudo make install -C build

If your computer has multiple cores, you may use the `make` command with the `-j` flag to use all your cores, or `-j` followed by a number to use the specified number of cores (for example, `-j4` to use 4 cores). You may perform the entire installation at once with:

    sudo cmake -B build && sudo make install -j -C build

Using lesser cores than you have is recommended so that your computer still has some headroom for other tasks. Thus, it may be a good idea for you to modify the above snippet, appending your preferred number of cores to `-j`.

### Troubleshooting

- If compilation suddenly stopped working for you, you may remove the `build` directory with `rm -rf build` and run the installation commands again.
- If running `cmake -B build` gives you an error, it is likely that a dependency was not installed correctly. Rerun the appropriate commands above to install the required dependencies and try again.
- If you get an error which mentions permission issues, try running the command that gave you the error with `sudo` prefixed, as done with `sudo make install -C build` above.
- If `cmake` complains that it could not find a particular library, please try uninstalling it and installing it again. After reinstalling it, please also unlink it and link it with the `--force` flag.
- If you get errors related to CGAL, it's likely you have another version of CGAL installed on your computer that is getting chosen instead of the one contained as a submodule within this repository. It's also possible that when cloning this repository, the `--recurse-submodule` flag was missing. Try running `git submodule init` and `git submodule update` in the root directory of the repository.

### Usage

Run the following command (replace `your-geojson-file.geojson` file with your geographic data and `your-csv-file.csv` with your visual variables file, containing target areas for each geographic region):

        cartogram your-geojson-file.geojson your-csv-file.csv

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

**You may find sample GeoJSON (containing geographic data) and CSV (containing information about target areas, colors and other visual variables) files in the `cartogram-cpp/sample_data` directory.**

### Testing

If you'd like to contribute to the project, please run our tests after you make any changes.

To run the unit tests, execute the following command:

    ctest --verbose

To learn more about the tests, you may go to the `cartogram-cpp/tests` directory and read the `README.md` file.

Additionally, you may go to the `cartogram-cpp/tests` directory and run the following command:

     bash stress_test.sh

### Uninstallation

Go to the `cartogram-cpp` directory in your preferred terminal and execute the following command:

    sudo make uninstall -C build

Upon successful uninstallation, the following will be outputted:

    > Built target uninstall

Further, running `cartogram` should no longer work.
