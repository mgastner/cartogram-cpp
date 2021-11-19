# lots of warnings
add_compile_options(-Wall -Wextra -pedantic)

# As per https://github.com/CGAL/cgal/wiki/How-to-use-CGAL-with-CMake-or-your-own-build-system
cmake_minimum_required(VERSION 3.1)
project(cartogram LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(THREADS_PREFER_PTHREAD_FLAG ON)
set(CMAKE_BUILD_TYPE Release)

set(CMAKE_CXX_COMPILER "g++-10")
set(CMAKE_C_COMPILER "gcc-10")

# Finding appropriate packages
find_package(CGAL 5.3 REQUIRED)
find_package(Boost 1.70 COMPONENTS program_options REQUIRED)
find_package(Threads REQUIRED)
find_package(PkgConfig REQUIRED)
pkg_search_module(FFTW REQUIRED fftw3 IMPORTED_TARGET)
include_directories(PkgConfig::FFTW)
link_libraries(PkgConfig::FFTW)

# Adding executables needed
add_executable(
  cartogram
  ../cartogram_generator/main.cpp
  ../cartogram_generator/geo_div.cpp
  ../cartogram_generator/cartogram_info.cpp
  ../cartogram_generator/inset_state.cpp
  ../cartogram_generator/read_csv.cpp
  ../cartogram_generator/rescale_map.cpp
  ../cartogram_generator/write_eps.cpp
  ../cartogram_generator/colors.cpp
  ../cartogram_generator/auto_color.cpp
  ../cartogram_generator/check_topology.cpp
  ../cartogram_generator/matrix.cpp
  ../cartogram_generator/densify.cpp
  ../cartogram_generator/densification_points.cpp
  ../cartogram_generator/read_geojson.cpp
  ../cartogram_generator/fill_with_density.cpp
  ../cartogram_generator/blur_density.cpp
  ../cartogram_generator/flatten_density.cpp
  ../cartogram_generator/ft_real_2d.cpp
  ../cartogram_generator/interpolate_bilinearly.cpp
  ../cartogram_generator/project.cpp
  ../cartogram_generator/write_geojson.cpp
  ../cartogram_generator/albers_projection.cpp
  ../cartogram_generator/simplify_inset.cpp)

# Linking appropriate libraries required
target_link_libraries(
  cartogram
  ${Boost_LIBRARIES}
  CGAL::CGAL
  ${CMAKE_THREAD_LIBS_INIT}
  PkgConfig::FFTW
)
