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

<a name="docker_on_mac"></a>
## Setting up a Docker Ubuntu Container (Mac instructions)

\* **This set up currently only allows for Vim as the text editor.**

1. If you have not already, download Docker Desktop from https://www.docker.com/products/docker-desktop.

2. Start Docker Desktop.

3. Change directories to your desired folder.

```
$ cd ~
```

4. Pull the Ubuntu Docker image.

```
$ docker pull ubuntu
```

5. View your Docker images.
```
$ docker image ls
```

6. Create a container named `ubuntu-cartogram_cpp` based on the Ubuntu image. The `-i` flag allows for an interactive container while the `-t` flag allows for terminal support within the container.
```
$ docker create -it --name=ubuntu-cartogram_cpp ubuntu
```

7. View your Docker containers and their relevant information.
```
$ docker ps -as
```

8. Start the `ubuntu-cartogram_cpp` container. You will enter the root directory of this Ubuntu container in an environment where you can run your normal terminal commands (e.g., `cd`, `cp`, `mv` `pwd`).
```
$ docker start -i ubuntu-cartogram_cpp
```

9. To start a second (or third, fourth, etc.) terminal window in this container, open a new terminal tab and run the following command.
```
$ docker exec -it ubuntu-cartogram_cpp bash
```

10. To exit the container, run the following command.
```
root@<number>:/# exit
```
11. Start the container and install the necessary dependencies for cartogram_cpp. You will need to use Vim as your text-editor (basic Vim instructions [here](https://opensource.com/article/19/3/getting-started-vim)).
```
$ docker start -i ubuntu-cartogram_cpp
root@<number>:/#
root@<number>:/# apt-get update
root@<number>:/# apt-get install sudo
root@<number>:/# apt-get install git
root@<number>:/# apt install vim
root@<number>:/# apt-get install cmake
root@<number>:/# apt install build-essential
root@<number>:/# apt-get install llvm
root@<number>:/# apt-get install libboost-all-dev
root@<number>:/# apt install g++-10
root@<number>:/# apt-get install libcgal-dev
root@<number>:/# apt-get install libomp-dev
```

12. Download and copy the nlohmann JSON parser library to your Ubuntu container.
    1. Go to https://github.com/nlohmann/json
    2. Click on "Code" -> "Download Zip"
    3. Go to your Downloads folder.
    4. Move the downloaded file from your Downloads folder to any other desired folder (not in the Ubuntu container) (e.g., root directory: `~`).
    5. Open a new separate terminal window and navigate to that desired folder.
    ```
    $ cd ~
    ```
    4. Copy the downloaded file to your Ubuntu container.
    ```
    $ docker cp json-develop/ ubuntu-cartogram_cpp:/home/json-develop/
    ```

13. To install nlohmann's JSON parser library, go back to your Ubuntu container terminal window and run the following commands.
```
root@<number>:/# cd home/json-develop/
root@<number>:/home/json-develop# cmake .
root@<number>:/home/json-develop# make
root@<number>:/home/json-develop# make install
```

14. Download and copy the FFTW library to your Ubuntu container.
    1. Go to [FFTW's website](http://www.fftw.org/download.html "FFTW Downloads Page").
    2. Install the latest version of FFTW (http)
    3. Go to Downloads folder
    4. Move the downloaded file from the Downloads folder to any desired folder (e.g., root directory: `~`).
    5. Open a new separate terminal window and navigate to that desired folder.
    ```
    $ cd ~
    ```
    4. Copy the downloaded file to your Ubuntu container.
    ```
    $ docker cp fftw-3.3.9.tar ubuntu-cartogram_cpp:/home/fftw-3.3.9.tar
    ```

15. To unarchive and install the fftw3 library, go back to your Ubuntu container terminal window and run the following commands.
```
root@<number>:/# cd /home/
root@<number>:/home# tar -xf fftw-3.3.9.tar
root@<number>:/home# rm fftw-3.3.9.tar
root@<number>:/home# cd fftw-3.3.9
root@<number>:/home/fftw-3.3.9# ./configure
root@<number>:/home/fftw-3.3.9# make
root@<number>:/home/fftw-3.3.9# sudo make install
```

16. In your Ubuntu container, follow the instructions [here](https://docs.github.com/en/github/getting-started-with-github/set-up-git) to set up Git and the Linux instructions [here](https://docs.github.com/en/github/authenticating-to-github/connecting-to-github-with-ssh) to connect to GitHub with SSH. Alternatively, you may use the following steps (though still referencing this GitHub guide).

    1. Set up Git. (Git was already installed in step 11)
    ```
    root@<number>:/# git config --global user.name "Your name here"
    root@<number>:/# git config --global user.email "your_email@example.com"
    ```

    2. Generate a new SSH key. Refer [here](https://docs.github.com/en/github/authenticating-to-github/generating-a-new-ssh-key-and-adding-it-to-the-ssh-agent) for how to fill in the fields when prompted.
    ```
    root@<number>:/# ssh-keygen -t ed25519 -C "your_email@example.com"
    ```
    3. Add your SSH key to the ssh-agent.
    ```
    root@<number>:/# eval "$(ssh-agent -s)"
    root@<number>:/# ssh-add ~/.ssh/id_ed25519
    ```

    4. Adding a new SSH key to your GitHub account.
        1. Open `.ssh/id_ed25519.pub` using Vim and manually copy its contents.
        ```
        root@<number>:/# vim ~/.ssh/id_ed25519.pub
        ```
        2. Exit Vim by typing `esc`, `:q` and enter.
        3. Create a new SSH key on GitHub and paste the contents into the 'key' section.

17. Clone the `cartogram_cpp` repository and begin developing in Vim as per normal.
```
root@<number>:/home# git clone git@github.com:mgastner/cartogram_cpp.git
root@<number>:/home# cd cartogram_cpp/build
root@<number>:/home/cartogram_cpp/build# cmake .
root@<number>:/home/cartogram_cpp/build# make
root@<number>:/home/cartogram_cpp/build# ./cartogram -h
```

## Installing Dependencies on macOS

**As of 2021-May-21, these instructions do not work on all Macs. We are trying to find the bug. In the meantime, please follow the [instructions for a Docker installation on a Mac](#docker_on_mac).**

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

`brew install libomp`

#### Installing LLVM Compilers

[LLVM Homepage](https://llvm.org)
`brew install llvm`

#### Installing FFTW3

[FFTW Homepage](http://www.fftw.org)

`brew install fftw`

#### Alternatively, install all, except nlhomann's JSON parser, at once:

`brew install fftw llvm libomp cgal`

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
