# Setting up CMake Project
cmake_minimum_required(VERSION 3.1)
project(cartogram LANGUAGES CXX)

# Setting compile options.
set(CMAKE_CXX_STANDARD 20)
set(THREADS_PREFER_PTHREAD_FLAG ON)
set(CMAKE_BUILD_TYPE Release)

# Lots of warnings.
add_compile_options(-Wall -Wextra -pedantic)

# Finding appropriate packages.
find_package(OpenMP REQUIRED)
find_package(CGAL 5.3 REQUIRED)
find_package(Boost 1.70 REQUIRED)
find_package(Threads REQUIRED)
find_package(PkgConfig REQUIRED)

# Finding and linking FFTW using PkgConfig.
pkg_search_module(FFTW REQUIRED fftw3 IMPORTED_TARGET)
include_directories(PkgConfig::FFTW)
link_libraries(PkgConfig::FFTW)

# Adding C++ files needed for compilation.
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
  ../cartogram_generator/simplify_inset.cpp
  ../cartogram_generator/write_geojson.cpp
  ../cartogram_generator/albers_projection.cpp
  ../cartogram_generator/parse_arguments.cpp)

# Linking appropriate libraries required.
if(APPLE)

  # Using GNU gcc-11, installed via homebrew, instaed of Apple clang.
  set(CMAKE_CXX_COMPILER "/usr/local/bin/g++-11")
  set(CMAKE_C_COMPILER "/usr/local/bin/gcc-11")

  # Adding relevant flag for OpenMP.
  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")

  # Linking appropriate libraries to cartogram target.
  target_link_libraries(
    cartogram
    ${Boost_LIBRARIES}
    CGAL::CGAL
    ${CMAKE_THREAD_LIBS_INIT}
    PkgConfig::FFTW
  )
elseif(UNIX)

  # Using GNU gcc-11.
  set(CMAKE_CXX_COMPILER "g++-11")
  set(CMAKE_C_COMPILER "gcc-11")

  # Linking appropriate libraries to cartogram target.
  target_link_libraries(
    cartogram
    ${Boost_LIBRARIES}
    CGAL::CGAL
    ${CMAKE_THREAD_LIBS_INIT}
    PkgConfig::FFTW
    OpenMP::OpenMP_CXX
  )
endif()

# Providing make with install target.
install(TARGETS cartogram DESTINATION bin)

# Providing make with uninstall target.
# TODO: Polish
# From https://stackoverflow.com/questions/41471620/cmake-support-make-uninstall
add_custom_target("uninstall")
add_custom_command(
    TARGET "uninstall"
    POST_BUILD
    COMMENT "Uninstalling cartogram..."
    COMMAND xargs rm -vf < install_manifest.txt || echo Nothing in
            install_manifest.txt to be uninstalled!
)
