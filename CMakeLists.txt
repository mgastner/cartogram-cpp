cmake_minimum_required(VERSION 3.27)
project(cartogram LANGUAGES CXX)

# ========== Project Setup ==========
set(CMAKE_CXX_STANDARD 20)
set(THREADS_PREFER_PTHREAD_FLAG ON)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# Default build type
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING "" FORCE)
endif()

# ========== Dependencies Setup ==========

# Boost
find_package(Boost REQUIRED)

# Matplot++
find_package(Matplot++ REQUIRED)

# PkgConfig, fftw, and cairo
find_package(PkgConfig REQUIRED)
pkg_search_module(fftw REQUIRED fftw3 IMPORTED_TARGET)
pkg_search_module(CAIRO REQUIRED CAIRO IMPORTED_TARGET)

# ========== Compiler Setup ==========
if(APPLE)
  set(LLVM_BASE_PATH "/usr/local/opt/llvm@17/bin/")

  if(${CMAKE_HOST_SYSTEM_PROCESSOR} STREQUAL "arm64")
    set(LLVM_BASE_PATH "/opt/homebrew/opt/llvm@17/bin/")
  endif()

  set(CMAKE_CXX_COMPILER "${LLVM_BASE_PATH}clang++")
  set(CMAKE_C_COMPILER "${LLVM_BASE_PATH}clang")
elseif(UNIX)
  set(CMAKE_CXX_COMPILER "g++-11")
  set(CMAKE_C_COMPILER "gcc-11")
endif()

# ========== Source Files ==========
file(GLOB_RECURSE CARTOGRAM_SOURCES "src/*.cpp")
add_executable(cartogram ${CARTOGRAM_SOURCES})

# ========== Include Directories ==========
if(Boost_FOUND)
  target_include_directories(cartogram PUBLIC ${Boost_INCLUDE_DIRS})
endif()

target_include_directories(cartogram PUBLIC
  ${PROJECT_SOURCE_DIR}/include
  ${Boost_INCLUDE_DIRS}
  PkgConfig::fftw
  PkgConfig::CAIRO
)

# ========== Compile Options ==========
if(CMAKE_COMPILER_IS_GNUCXX OR CMAKE_CXX_COMPILER_ID MATCHES "Clang")
  target_compile_options(cartogram PRIVATE -isystem ${Boost_INCLUDE_DIRS})
  target_compile_options(cartogram PRIVATE -ffp-contract=off)
elseif(MSVC)
  target_compile_options(cartogram PRIVATE /external:I ${Boost_INCLUDE_DIRS})
endif()

# Compiler warnings
target_compile_options(cartogram PRIVATE -Wall -Wextra -pedantic -Wno-deprecated-declarations)

# ========== Linking Libraries ==========
target_link_libraries(cartogram
  PkgConfig::fftw
  PkgConfig::CAIRO
  Matplot++::matplot
)

# ========== Installation ==========
install(TARGETS cartogram DESTINATION bin)

# Uninstall target
add_custom_target("uninstall")
add_custom_command(
  TARGET "uninstall"
  POST_BUILD
  COMMENT "Uninstalling cartogram..."
  COMMAND xargs rm -vf < install_manifest.txt || echo "Nothing in install_manifest.txt to be uninstalled!"
)
