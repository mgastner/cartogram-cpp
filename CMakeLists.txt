cmake_minimum_required(VERSION 3.22)
project(cartogram LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# ----------------------------- Third–party -----------------------------
find_package(Boost REQUIRED CONFIG)
find_package(CGAL REQUIRED CONFIG)
find_package(nlohmann_json REQUIRED CONFIG)
find_package(argparse REQUIRED CONFIG)
find_package(vincentlaucsb-csv-parser REQUIRED CONFIG)
find_package(indicators REQUIRED CONFIG)
find_package(PkgConfig REQUIRED)
pkg_search_module(FFTW REQUIRED IMPORTED_TARGET fftw3)
pkg_search_module(CAIRO REQUIRED IMPORTED_TARGET cairo)

# ----------------------------- Main exe  ------------------------------
file(GLOB_RECURSE CARTOGRAM_SOURCES CONFIGURE_DEPENDS src/*.cpp)
add_executable(cartogram ${CARTOGRAM_SOURCES})

target_link_libraries(cartogram
  PRIVATE
  CGAL::CGAL
  Boost::boost
  nlohmann_json::nlohmann_json
  argparse::argparse
  vincentlaucsb-csv-parser::vincentlaucsb-csv-parser
  indicators::indicators
  PkgConfig::FFTW
  PkgConfig::CAIRO
)

target_compile_options(cartogram PRIVATE
  -ffp-contract=off
  -Wall -Wextra -Wpedantic -Wshadow -Wold-style-cast -Wunreachable-code
)

target_include_directories(cartogram
  PUBLIC ${PROJECT_SOURCE_DIR}/include
)

install(TARGETS cartogram RUNTIME DESTINATION bin)
