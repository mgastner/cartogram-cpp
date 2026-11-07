cmake_minimum_required(VERSION 3.30)
project(cartogram LANGUAGES CXX)

enable_testing()

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
find_package(jsoncons REQUIRED)

pkg_search_module(FFTW REQUIRED IMPORTED_TARGET fftw3)

# ----------------------------- Compiler flags  ------------------------------
# Shared by cartogram_lib and conformality_core
set(CARTOGRAM_STRICT_FLAGS
  -ffp-contract=off
  -Werror
  -Wall -Wextra -Wpedantic -Wshadow -Wold-style-cast -Wunreachable-code
  -pedantic-errors -Wformat=2 -Wformat-security -Wcast-qual -Wcast-align -Wconversion
  -Wsign-conversion -Wfloat-equal -Woverloaded-virtual -Wno-unknown-pragmas
  -Wmissing-declarations -Werror=return-type -Wnon-virtual-dtor
)

# ----------------------------- Cartogram Library  ---------------------------
file(GLOB_RECURSE CARTOGRAM_SOURCES CONFIGURE_DEPENDS src/*.cpp)

# Exclude main.cpp from the list of sources
list(FILTER CARTOGRAM_SOURCES EXCLUDE REGEX ".*/main\\.cpp$")

# The conformality core is built as its own library (see below)
list(FILTER CARTOGRAM_SOURCES EXCLUDE REGEX ".*/src/conformality/.*")

add_library(cartogram_lib STATIC ${CARTOGRAM_SOURCES})
target_include_directories(cartogram_lib
  PUBLIC ${PROJECT_SOURCE_DIR}/include
)

target_link_libraries(cartogram_lib
  PUBLIC
  CGAL::CGAL
  boost::boost
  nlohmann_json::nlohmann_json
  argparse::argparse
  vincentlaucsb-csv-parser::vincentlaucsb-csv-parser
  indicators::indicators
  PkgConfig::FFTW
  jsoncons
)

if(CMAKE_BUILD_TYPE STREQUAL Release)
  target_compile_options(cartogram_lib PUBLIC -O3)
else()
  target_compile_options(cartogram_lib PUBLIC -O1)

  # Ignore CGAL preconditions in Debug builds (false positives)
  target_compile_definitions(cartogram_lib
    PUBLIC CGAL_NO_PRECONDITIONS
  )
endif()

target_compile_options(cartogram_lib PUBLIC ${CARTOGRAM_STRICT_FLAGS})

if(CMAKE_BUILD_TYPE STREQUAL Release)
  include(CheckCXXCompilerFlag)
  check_cxx_compiler_flag("-march=native" HAS_MARCH_NATIVE)
  check_cxx_compiler_flag("-mtune=native" HAS_MTUNE_NATIVE)

  if(HAS_MARCH_NATIVE)
    target_compile_options(cartogram_lib PUBLIC -march=native)
  elseif(HAS_MTUNE_NATIVE)
    target_compile_options(cartogram_lib PUBLIC -mtune=native)
  endif()

  include(CheckIPOSupported)
  check_ipo_supported(RESULT ipo_ok OUTPUT ipo_msg)

  if(CMAKE_CXX_COMPILER_ID STREQUAL "Clang" AND CMAKE_SYSTEM_NAME STREQUAL "Linux")
    if(ipo_ok)
      find_program(LLD_PATH NAMES lld ld.lld)

      if(NOT LLD_PATH)
        message(WARNING "lld not found; disabling IPO. Install with: sudo apt install lld")
        set(ipo_ok FALSE)
      endif()
    endif()
  endif()

  if(ipo_ok)
    set_property(TARGET cartogram_lib PROPERTY INTERPROCEDURAL_OPTIMIZATION ON)
    target_link_options(cartogram_lib PUBLIC -flto=full)

    if(CMAKE_CXX_COMPILER_ID STREQUAL "Clang" AND CMAKE_SYSTEM_NAME STREQUAL "Linux")
      message(STATUS "Using LLD with full LTO + whole-program vtables: ${LLD_PATH}")
      target_link_options(cartogram_lib PUBLIC -fuse-ld=lld)
      target_compile_options(cartogram_lib
        PUBLIC
        -flto=full
        -fwhole-program-vtables
        -fvirtual-function-elimination)
    endif()
  else()
    message(STATUS "IPO/LTO not available: ${ipo_msg}")
  endif()
endif()

target_include_directories(cartogram_lib
  PUBLIC ${PROJECT_SOURCE_DIR}/include
)

# ----------------------------- Conformality core  ---------------------------
# Numerical core of the Jacobian-preserving conformality postprocessor. It
# has no dependency on InsetState, CGAL, or the other third-party libraries.
file(GLOB_RECURSE CONFORMALITY_SOURCES CONFIGURE_DEPENDS
  src/conformality/*.cpp)

add_library(conformality_core STATIC ${CONFORMALITY_SOURCES})
target_include_directories(conformality_core
  PUBLIC ${PROJECT_SOURCE_DIR}/include
)

if(CMAKE_BUILD_TYPE STREQUAL Release)
  target_compile_options(conformality_core PUBLIC -O3)
  if(HAS_MARCH_NATIVE)
    target_compile_options(conformality_core PUBLIC -march=native)
  elseif(HAS_MTUNE_NATIVE)
    target_compile_options(conformality_core PUBLIC -mtune=native)
  endif()
else()
  target_compile_options(conformality_core PUBLIC -O1)
endif()

target_compile_options(conformality_core PUBLIC ${CARTOGRAM_STRICT_FLAGS})

# Optional multithreading of the core's node loops. Results do not depend on
# the thread count (see include/conformality/parallel.hpp). OpenMP is linked
# privately, so the rest of cartogram_lib is still compiled without it. With
# Apple Clang, pass -DOpenMP_ROOT=$(brew --prefix libomp).
option(CONFORMALITY_OPENMP "Parallelize the conformality core with OpenMP" OFF)
if(CONFORMALITY_OPENMP)
  find_package(OpenMP REQUIRED COMPONENTS CXX)
  target_link_libraries(conformality_core PRIVATE OpenMP::OpenMP_CXX)
endif()

# The cartogram library hosts the adapter that feeds the core
target_link_libraries(cartogram_lib PUBLIC conformality_core)

# ----------------------------- Main executable  -----------------------------
add_executable(cartogram src/main.cpp)
target_link_libraries(cartogram
  PRIVATE cartogram_lib
)
install(TARGETS cartogram RUNTIME DESTINATION bin)

# Independent area check of a logged conformality correction (research tool)
add_executable(conformality_check research/tools/conformality_check.cpp)
target_link_libraries(conformality_check
  PRIVATE conformality_core nlohmann_json::nlohmann_json
)

# ----------------------------- Tests  -----------------------------
add_subdirectory(tests)
