cmake_minimum_required(VERSION 3.30)

project(cartogram LANGUAGES CXX)

# ========== Project Setup ==========
set(CMAKE_CXX_STANDARD 20)
set(THREADS_PREFER_PTHREAD_FLAG ON)
set(CMAKE_COLOR_DIAGNOSTICS ON)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# Assume development build by default
set(RELEASE_TAG "development" CACHE STRING "Release tag for the build")

# Default build type
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING "" FORCE)
endif()

# ========== Dependencies Setup ==========
find_package(Boost REQUIRED CONFIG)
find_package(cgal REQUIRED CONFIG)
find_package(nlohmann_json REQUIRED CONFIG)
find_package(PkgConfig REQUIRED)
pkg_search_module(FFTW REQUIRED IMPORTED_TARGET fftw3)

# ========== Source Files ==========
file(GLOB_RECURSE CARTOGRAM_SOURCES "src/*.cpp")
add_executable(cartogram ${CARTOGRAM_SOURCES})

target_compile_definitions(cartogram PRIVATE RELEASE_TAG="${RELEASE_TAG}")

# ========== Include Directories ==========
target_include_directories(cartogram
  PUBLIC
  ${PROJECT_SOURCE_DIR}/include
)

target_include_directories(cartogram
  SYSTEM PUBLIC
  ${PROJECT_SOURCE_DIR}/external
)

# ========== Compile Options ==========
target_compile_options(cartogram PRIVATE -ffp-contract=off -O3 -Wno-unknown-pragmas)

# Compiler warnings
target_compile_options(cartogram PRIVATE
  -Wall # Enable all warnings
  -Wextra # Enable extra warnings
  -Wpedantic # Enable pedantic warnings
)

# ========== Linking Libraries ==========
target_link_libraries(cartogram
  CGAL::CGAL
  boost::boost
  nlohmann_json::nlohmann_json
  PkgConfig::FFTW
)

# ========== Installation ==========
install(TARGETS cartogram DESTINATION /usr/local/bin)

# Uninstall target
add_custom_target("uninstall")
add_custom_command(
  TARGET "uninstall"
  POST_BUILD
  COMMENT "Uninstalling cartogram..."
  COMMAND xargs rm -vf < install_manifest.txt || echo "Nothing in install_manifest.txt to be uninstalled!"
)
