# cartogram_cpp
Cartogram generator in C++

## Installing Dependencies on Ubuntu

#### Installing nlohmann's JSON parser
1. Go to https://github.com/nlohmann/json
2. Click on "Code" -> "Download Zip"
3. Go to Downloads folder
4. Unzip the file you just installed (you can use the `unzip` command)
5. Go into the newly created unzipped folder json-develop (you can use the `cd` command)

```
cmake .
make
sudo make install
```

<!-- #### Installing CGAL

[CGAL Homepage](https://www.cgal.org/)

`sudo apt-get install libcgal-dev` -->

#### Installing CGAL

[CGAL Homepage](https://www.cgal.org/)

CGAL Version 5.3 is currently unavailable from apt.
Please follow the instructions on the CGAL website to build from source.

You may download the latest release [here](https://github.com/CGAL/cgal/releases).
You may find installation instructions [here](https://doc.cgal.org/latest/Manual/usage.html#title4).


For posterity: Once version 5.3 is available through apt-get (you may check [here](https://packages.ubuntu.com/search?keywords=libcgal-dev&searchon=names&suite=impish&section=all)), you may run the following command to install it.

```
sudo apt-get install libcgal-dev
```

#### Installing OpenMP

[OpenMP Homepage](https://www.openmp.org/)

```
sudo apt-get install libomp-dev
```

#### Installing FFTW3
1. Go to [FFTW's website](http://www.fftw.org/download.html "FFTW Downloads Page").
2. Install the latest version of FFTW
3. Unarchive the file with: `tar zxvf fftw-3.3.9.tar.gz` (Note: the version number may be different).
4. Go to the directory with: `cd fftw-3.3.9`
```
./configure
make
7sudo make install
```

## Setting up dependencies on macOS (ARM & x86)

### Intel-Only Instructions (x86, Macs released before 2020)

1. Install [homebrew](brew.sh) using:
```
/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
```

2. Install gcc-10, icu4c, pkg-config, wget, boost, fftw, cgal, nlohmann-json and cmake.
```
brew install gcc@11 icu4c pkg-config wget boost fftw cgal nlohmann-json cmake
```

---

### ARM-Only Instructions (M1, M1 Pro, M1 Max, etc.)

1. Go to your applications folder.
2. Right-click on your Terminal and duplicate it.
3. Rename your newly duplicated Terminal to `x86 Terminal` (or something else of your choice).
4. Right-click on the new Terminal, select Get Info.
5. On this page, make sure to select `Open using Rosetta`.
6. You now have an x86 terminal, which you can use to build other x86 binaries as well.
7. Install [homebrew](brew.sh) using:
```
/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
```
8. Open your `zshrc` using:
```
touch ~/.zshrc
open -a TextEdit ~/.zshrc
```
9. Insert the following line:
```
alias brew86='arch -x86_64 /usr/local/Homebrew/bin/brew'
```
and save the file.

10. Make sure Zsh knows you've updated your `zshrc` with:
```
source ~/.zshrc
```
You may confirm that you followed the instructions correctly by executing.
```
type brew86
```

The output for the above command should be:
> brew86 is an alias for arch -x86_64 /usr/local/Homebrew/bin/brew

10. Finally, install the required dependencies by running:
```
brew86 install gcc@11 icu4c pkg-config wget boost fftw cgal nlohmann-json cmake
```

11. You must use the `x86 Terminal` that you created in steps 1-6 to compile and run the program.

## Compilation and Usage (Ubuntu or macOS)

#### Compilation

Go to the `cartogram_cpp` directory and execute the following commands.

```
cd build
cmake .
make
```

#### Usage

1. Go to the build directory using `cd build`.
2. Replace `your-geojson-file.geojson` file with your geographic data and `your-csv-file.csv` with your visual variables file, and run the following command:
```
./cartogram your-geojson-file.geojson -v your-csv-file.csv`
```

- The first positional argument's input is a GeoJSON or JSON file, in the standard GeoJSON format.
- The `-v` flag accepts a `.csv` file with data about target areas.

*Note: use the `-h` flag to display more options*

The csv file should be in the following format:

| NAME_1        | Data (eg: Population)| Color   |
| :------------ |:---------------------| :-------|
| Bruxelles     | 1208542              | #e74c3c |
| Vlaanderen    | 6589069              | #f1c40f |
| Wallonie      | 3633795              | #34495e |

- `NAME_1` should be the same as the identifying property's name in the GeoJSON. The rows should also have the same data as is present in the identifying property.
- `Data` contains the data you would like your cartogram to based on.
- `Color` is the color you would like the geographic region to be.

Colors may be represented in the following manner:
1. `cornflowerblue`: html color codes supported by `CSS3` (case-insensitive), full list of supported colors may be found in the "Extended colors" section of [web colors](https://en.wikipedia.org/wiki/Web_colors).
2. `"rgb(255, 0, 120)"` or `rgb(255 0 120)` or `"255, 0, 120"` or `255 0 120`: red, green and blue values out of 255
3. `#e74c3c`: hex code of color, must start with `#`
