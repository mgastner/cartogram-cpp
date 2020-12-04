cmake_minimum_required(VERSION 2.8.2)

project(findFFTW-download NONE)

include(ExternalProject)

ExternalProject_Add(findFFTW_download
    GIT_REPOSITORY    "https://github.com/egpbos/findfftw.git"
    CONFIGURE_COMMAND ""
    BUILD_COMMAND     ""
    INSTALL_COMMAND   ""
    TEST_COMMAND      ""
    SOURCE_DIR        "/mnt/d/projects/github/cartogram_cpp_master/cartogram_generator/findFFTW-src"
    BINARY_DIR        ""
    INSTALL_DIR       ""
)
