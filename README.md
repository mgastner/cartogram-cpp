# cartogram_cpp
Cartogram generator in C++

## Installing Dependencies on Ubuntu

#### Installing nlohmann's JSON parser
1. Go to https://github.com/nlohmann/json
2. Click on "Code" -> "Download Zip"
3. Go to Downloads folder
4. Unzip the file you just installed (you can use the `unzip` command)
5. Go into the newly created unzipped folder json-develop (you can use the `cd` command)
6. `cmake .`
7. `make`
8. `sudo make install`

#### Installing CGAL

[CGAL Homepage](https://www.cgal.org/)

`sudo apt-get install libcgal-dev`

#### Installing OpenMP

[OpenMP Homepage](https://www.openmp.org/)

`sudo apt-get install libomp-dev`


#### Installing FFTW3
1. Go to [FFTW's website](http://www.fftw.org/download.html "FFTW Downloads Page").
2. Install the latest version of FFTW
3. Unarchive the file with: `tar zxvf fftw-3.3.9.tar.gz`
4. Go to the directory with: `cd fftw-3.3.9`
5. `./configure`
6. `make`
7. `sudo make install`

## Setting up dependencies on macOS (ARM & x86)

### ARM Only Instructions (M1, M1 Pro, M1 Max, etc.)

1. Go to your applications folder.
2. Right-click on your Terminal and duplicate it.
3. Rename your newly duplicated Terminal to `x86 Terminal` (or something else of your choice).
4. Right-click on the new Terminal, select Get Info.
5. On this page, make sure to select `Open using Rosetta`.
6. You now have a x86 terminal, which you can use to build other x86 binaries as well.
7. Install [homebrew](brew.sh) using:
```
/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
```
8. Open your zshrc using:
```
touch ~/.zshrc
open -a TextEdit ~/.zshrc
```
9. Insert the following line:
```
alias brew86='arch -x86_64 /usr/local/Homebrew/bin/brew'
```
and save the file.
10. Make sure your Zsh knows you've updated your zshrc with:
```
source ~/.zshrc
```
You can confirm executing the above command and then trying `brew86`.
11. You may now start from step 2 of the General Instructions.

##### **Please make sure you replace `brew` with `brew86` in the steps below**. Otherwise, you will NOT be able to build the binary on your ARM machine.

Further, make sure you *skip* step 1, as you would have already installed homebrew at this point.

### General Instructions (If you have an Intel processor, start here)

1. Install [homebrew](brew.sh) using:
```
/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
```

2. Install gcc-10, icu4c, pkg-config and wget.

`brew install gcc@10 icu4c pkg-config wget`

4. Build Boost from Source using gcc-10

Please ensure you do not have an existing installation of boost. If you usually use homebrew to install packages, you can ensure that you uninstall it with: `brew uninstall boost`.

Run the following commands, one after another, to build boost from source:
```
wget https://boostorg.jfrog.io/artifactory/main/release/1.75.0/source/boost_1_75_0.tar.gz
tar -xf boost_1_75_0.tar.gz
cd boost_1_75_0
CC=gcc-10 CXX=g++-10 ./bootstrap.sh --with-icu=/usr/local/opt/icu4c --without-libraries=python,mpi
./b2 headers
echo "using darwin : : g++-10 ;\n" > user-config.jam
./b2 -d2 -j8 --layout=tagged-1.66 --user-config=user-config.jam -sNO_LZMA=1 -sNO_ZSTD=1 install threading=multi,single
```
*Note: This may take up to 15 minutes, depending on your computer.*

Thanks to [this page](https://githubmemory.com/repo/Homebrew/homebrew-core/issues/76645).

4. Install FFTW with our settings (without MPI or OpenMP)

[FFTW Homepage](http://www.fftw.org)

Navigate to the directory with our modified `fftw.rb` file `build/fftw/fftw.rb`, based on the original homebrew ruby script on homebrew-core [here](https://github.com/Homebrew/homebrew-core/blob/master/Formula/fftw.rb). Then, run the following command.
```
brew install --build-from-source fftw.rb
```

2. Install nlohmann's JSON parser and CGAL.

[nlohmann's GitHub Page](https://github.com/nlohmann/json)
[CGAL Homepage](https://www.cgal.org/)
```
brew tap nlohmann/json
brew install nlohmann-json cgal
```

## Compiling and running (Ubuntu or macOS)

#### Compile by running:

1. `cd ./build`
2. `cmake .`
3. `make`

#### To use the cartogram generator:

1. `cd ./build`
2. `./cartogram -g your-geojson-file.geojson -v your-csv-file.csv`

- The `-g` flag accepts a GeoJSON or JSON file, in the standard GeoJSON format.
- The `-v` flag accepts a .csv file with your target areas data.

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
