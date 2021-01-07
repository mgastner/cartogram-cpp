# cartogram_cpp
Cartogram generator in C++

## Note: This is still a work-in-progress and, hence, cartograms are not outputted at the moment.

## Installing Dependencies on Ubuntu

#### Installing nlohmann's JSON parser
1. Go to https://github.com/nlohmann/json
2. Click on "Code" -> "Download Zip"
3. Go to Downloads folder
4. `cmake .`
5. `make`
6. `sudo make install`

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

## Installing Dependencies on macOS

**Please ensure you have `Homebrew` installed. You may find instruction on [Homebrew's website](https://brew.sh "Homebrew Home Page") for the same.**

#### Installing nlohmann's JSON parser

[nlohmann's GitHub Page](https://github.com/nlohmann/json)

```
brew tap nlohmann/json
brew install nlohmann-json
```


#### Installing CGAL

[CGAL Homepage](https://www.cgal.org/)

`brew install cgal`

#### Installing OpenMP

[OpenMP Homepage](https://www.openmp.org/)

`brew install llvm libomp`


#### Installing FFTW3

[FFTW Homepage](http://www.fftw.org)

`brew install fftw`

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
