cmake_minimum_required(VERSION 3.30)
project(cartogram LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# ----------------------------- Third–party -----------------------------
find_package(Boost REQUIRED CONFIG)
find_package(cgal REQUIRED CONFIG)
find_package(nlohmann_json REQUIRED CONFIG)
find_package(argparse REQUIRED CONFIG)
find_package(vincentlaucsb-csv-parser REQUIRED CONFIG)
find_package(indicators REQUIRED CONFIG)
find_package(PkgConfig REQUIRED)
pkg_search_module(FFTW REQUIRED IMPORTED_TARGET fftw3)

# ----------------------------- Main exe  ------------------------------
file(GLOB_RECURSE CARTOGRAM_SOURCES CONFIGURE_DEPENDS src/*.cpp)
add_executable(cartogram ${CARTOGRAM_SOURCES})

target_link_libraries(cartogram
  PRIVATE
  CGAL::CGAL
  boost::boost
  nlohmann_json::nlohmann_json
  argparse::argparse
  vincentlaucsb-csv-parser::vincentlaucsb-csv-parser
  indicators::indicators
  PkgConfig::FFTW
)

target_compile_options(cartogram PRIVATE
  -O3
  -ffp-contract=off
  -Wall -Wextra -Wpedantic -Wshadow -Wold-style-cast -Wunreachable-code
)

include(CheckCXXCompilerFlag)
check_cxx_compiler_flag("-mcpu=native" HAS_MCPU_NATIVE)
check_cxx_compiler_flag("-march=native" HAS_MARCH_NATIVE)

if (HAS_MCPU_NATIVE)
    target_compile_options(cartogram PRIVATE -mcpu=native)
elseif (HAS_MARCH_NATIVE)
    target_compile_options(cartogram PRIVATE -march=native)
endif()


include(CheckIPOSupported)
check_ipo_supported(RESULT ipo_ok OUTPUT ipo_msg)

if(ipo_ok)
  set_property(TARGET cartogram PROPERTY INTERPROCEDURAL_OPTIMIZATION ON)
  target_link_options(cartogram PRIVATE -flto)
else()
  message(STATUS "IPO/LTO not available: ${ipo_msg}")
endif()

target_include_directories(cartogram
  PUBLIC ${PROJECT_SOURCE_DIR}/include
)

install(TARGETS cartogram RUNTIME DESTINATION bin)
