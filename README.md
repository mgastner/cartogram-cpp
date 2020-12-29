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

## Compiling and running on Ubuntu

#### Compile by running:

1. `cd ./build`
2. `cmake .`
3. `make`

#### To use the cartogram generator:

1. `cd ./build`
2. `./cartogram -g your-geojson-file.geojson -v your-csv-file.csv`


- The `-g` flag accepts a GeoJSON or JSON file, in the standard GeoJSON format.
- The `-v` flag accepts a .csv file with your target areas data.

The csv file should be in the following format:

| NAME_1        | Data (eg: Population)| Color   |
| :------------ |:---------------------| :-------|
| Bruxelles     | 1208542              | #e74c3c |
| Vlaanderen    | 6589069              | #f1c40f |
| Wallonie      | 3633795              | #34495e |

- `NAME_1` should be the same as the identifying property's name in the GeoJSON. The rows should also have the same data as is present in the identifying property.
- `Data` contains the data you would like your cartogram to based on.
- `Color` is the hex color code you would like the geographic region to be colored.
