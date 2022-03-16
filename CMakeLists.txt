# Setting up CMake Project
cmake_minimum_required(VERSION 3.1)
project(cartogram LANGUAGES CXX)

# Setting compile options.
set(CMAKE_CXX_STANDARD 20)
set(THREADS_PREFER_PTHREAD_FLAG ON)

# Setting CMake build options, adding path for binary.
set(CMAKE_BUILD_TYPE Release)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# Lots of warnings.
add_compile_options(-Wall -Wextra -pedantic)

# Finding appropriate packages.
find_package(OpenMP REQUIRED)
find_package(PkgConfig REQUIRED)

# Finding and linking FFTW using PkgConfig.
pkg_search_module(FFTW REQUIRED fftw3 IMPORTED_TARGET)
pkg_search_module(Cairo REQUIRED cairo IMPORTED_TARGET)
include_directories(PkgConfig::FFTW PkgConfig::Cairo)
link_libraries(PkgConfig::FFTW PkgConfig::Cairo)

# Adding C++ files needed for compilation.
add_executable(
  cartogram
  src/main.cpp
  src/geo_div.cpp
  src/cartogram_info.cpp
  src/read_csv.cpp
  src/rescale_map.cpp
  src/write_cairo.cpp
  src/write_eps.cpp
  src/colors.cpp
  src/check_topology.cpp
  src/matrix.cpp
  src/read_geojson.cpp
  src/intersection.cpp
  src/blur_density.cpp
  src/flatten_density.cpp
  src/ft_real_2d.cpp
  src/interpolate_bilinearly.cpp
  src/project.cpp
  src/simplify_inset.cpp
  src/write_geojson.cpp
  src/albers_projection.cpp
  src/parse_arguments.cpp
  src/inset_state/auto_color.cpp
  src/inset_state/fill_with_density.cpp
  src/inset_state/inset_state.cpp
  src/inset_state/densify.cpp
  src/inset_state/scanline_graph.cpp
  src/smyth_projection.cpp
  src/round_point.cpp
)

# Linking appropriate libraries required.
if(APPLE)

  # Using LLVM clang, installed via homebrew, instaed of Apple clang.
  if(${CMAKE_HOST_SYSTEM_PROCESSOR} STREQUAL "arm64")

    # Homebrew install path on arm64 (apple silicon, M1 etc.) machines.
    set(CMAKE_CXX_COMPILER "/opt/homebrew/opt/llvm/bin/clang++")
    set(CMAKE_C_COMPILER "/opt/homebrew/opt/llvm/bin/clang")
  else()

    # Homebrew install path on x86 machines (Intel).
    set(CMAKE_C_COMPILER "/usr/local/opt/llvm/bin/clang")
    set(CMAKE_CXX_COMPILER "/usr/local/opt/llvm/bin/clang++")
  endif()
elseif(UNIX)

  # Using GNU gcc-11.
  set(CMAKE_CXX_COMPILER "g++-11")
  set(CMAKE_C_COMPILER "gcc-11")
endif()

# Linking appropriate libraries to cartogram target.
target_link_libraries(
  cartogram
  PkgConfig::FFTW
  PkgConfig::Cairo
  OpenMP::OpenMP_CXX
)

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
