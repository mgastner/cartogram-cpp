# cartogram-cpp: Cartogram generator in C++ [![DOI](https://zenodo.org/badge/281575635.svg)](https://zenodo.org/badge/latestdoi/281575635) [![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](./LICENSE)

<p align="center">
<img src ="img/gocart_logo.svg" width="50%">
</p>

This program uses the fast flow-based method developed by Michael T. Gastner, Vivien Seguy, and Pratyush More. For more information, you may refer to the following [paper](https://www.pnas.org/content/115/10/E2156):

Gastner MT, Seguy V, More P. _Fast flow-based algorithm for creating density-equalizing map projections_. Proc Natl Acad Sci USA 115(10):E2156–E2164 (2018). <https://doi.org/10.1073/pnas.0400280101>

Data produced by code in this repository are subject to the MIT license found [here](./LICENSE) and should cite the aforementioned paper by Gastner et al. (2018).

## Dependencies

### macOS

#### Installing Homebrew

Install [homebrew](brew.sh) by running the following command:

    /bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"

#### Installing dependencies through Homebrew

Install llvm, pkg-config, boost, fftw, cgal, nlohmann-json, and cmake by running the following command:

    brew install llvm libomp pkg-config boost fftw cgal nlohmann-json cmake cairo

### Debian-based distributions (Ubuntu, Arch Linux etc.)

#### Installing GNU gcc-11

GNU gcc-11 is currently unavailable from apt by default. You may find installation instructions [here](https://lindevs.com/install-gcc-on-ubuntu/). Alternatively, you may run the following commands to install it:

    sudo apt install build-essential manpages-dev software-properties-common
    sudo add-apt-repository ppa:ubuntu-toolchain-r/test
    sudo apt update && sudo apt install gcc-11 g++-11

#### Installing nlohmann's JSON parser

1.  Go to <https://github.com/nlohmann/json>
2.  Click on "Code" -> "Download Zip".
3.  Go to Downloads folder.
4.  Unzip the file you just installed (you can use the `unzip` command).
5.  Go into the newly created unzipped folder json-develop (you can use the `cd` command).
6.  Run the following commands (you may copy and paste all of them at once):

        cmake .
        make
        sudo make install

#### Installing CGAL

[CGAL Homepage](https://www.cgal.org/)

CGAL Version 5.3 is currently unavailable from apt. Please follow the instructions on the CGAL website to build from source.

You may download the latest release [here](https://github.com/CGAL/cgal/releases). You may find installation instructions [here](https://doc.cgal.org/latest/Manual/usage.html#title4).

For posterity: Once version 5.3 is available through apt (you may check [here](https://packages.ubuntu.com/search?keywords=libcgal-dev&searchon=names&suite=impish§ion=all)), you may run the following command to install it.

    sudo apt install libcgal-dev

#### Installing OpenMP

[OpenMP Homepage](https://www.openmp.org/)

    sudo apt install libomp-dev

#### Installing FFTW3

1.  Go to [FFTW's website](http://www.fftw.org/download.html "FFTW Downloads Page").
2.  Install the latest version of FFTW.
3.  Unarchive the file with: `tar zxvf fftw-3.3.10.tar.gz` (Note: the version number may be different).
4.  Go to the directory with: `cd fftw-3.3.10`.
5.  Run the following commands (you may copy and paste all of them at once):

        ./configure
        make
        sudo make install

#### Installing CairoGraphics

[CairoGraphics Homepage](https://www.cairographics.org/)

    sudo apt install libcairo2-dev        


### Installation

Go to the `cartogram_cpp` directory in your preferred terminal and execute the following commands.

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

### Usage

Run the following command (replace `your-geojson-file.geojson` file with your geographic data and `your-csv-file.csv` with your visual variables file, containing target areas for each geographic region):

        cartogram your-geojson-file.geojson your-csv-file.csv

-   The first argument's input is a GeoJSON or JSON file, in the standard GeoJSON format.
-   The second argument's input is a `.csv` file with data about target areas.

_Note: use the `-h` flag to display more options._

The CSV file should be in the following format:

| NAME_1     | Data (e.g., Population) | Color   |
|:-----------|:------------------------|:--------|
| Bruxelles  | 1208542                 | #e74c3c |
| Vlaanderen | 6589069                 | #f1c40f |
| Wallonie   | 3633795                 | #34495e |

-   `NAME_1` should be the same as the identifying property's name in the GeoJSON. The rows should also have the same data as is present in the identifying property.
-   `Data` contains the data you would like your cartogram to based on.
-   `Color` is the color you would like the geographic region to be. Colors may be represented in the following manner:

    1.  `cornflowerblue`: html color codes supported by `CSS3` (case-insensitive), full list of supported colors may be found in the "Extended colors" section of [web colors](https://en.wikipedia.org/wiki/Web_colors).
    2.  `"rgb(255, 0, 120)"` or `rgb(255 0 120)` or `"255, 0, 120"` or `255 0 120`: red, green and blue values out of 255.
    3.  `#e74c3c`: hex code of color, must start with `#`.

**You may find sample GeoJSON (containing geographic data) and CSV (containing information about target areas, colors and other visual variables) files in the `cartogram_cpp/sample_data` directory.**

### Testing

If you'd like to contribute to the project, please run our test battery after you make any changes. You may do so by going to the `cartogram_cpp/tests` directory and running the following command:

        bash stress_test.sh

### Uninstallation

Go to the `cartogram_cpp` directory in your preferred terminal and execute the following command:

    sudo make uninstall -C build

Upon successful uninstallation, the following will be outputted:

    > Built target uninstall

Further, running `cartogram` should no longer work.
