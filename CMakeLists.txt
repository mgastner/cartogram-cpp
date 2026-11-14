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

# ----------------------------- Cartogram Library  ---------------------------
file(GLOB_RECURSE CARTOGRAM_SOURCES CONFIGURE_DEPENDS src/*.cpp)

# Exclude main.cpp from the list of sources
list(FILTER CARTOGRAM_SOURCES EXCLUDE REGEX ".*/main\\.cpp$")

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

target_compile_options(cartogram_lib PUBLIC
  -ffp-contract=off
  -Werror
  -Wall -Wextra -Wpedantic -Wshadow -Wold-style-cast -Wunreachable-code
  -pedantic-errors -Wformat=2 -Wformat-security -Wcast-qual -Wcast-align -Wconversion
  -Wsign-conversion -Wfloat-equal -Woverloaded-virtual -Wno-unknown-pragmas
  -Wmissing-declarations -Werror=return-type -Wnon-virtual-dtor
)

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
      # With -fuse-ld=lld, clang runs the ld.lld in its own directory if
      # there is one, and otherwise the one on the PATH. Search in the same
      # order, so that LTO is not disabled when only a versioned package
      # such as lld-22 is installed, and the version found is the one used.
      file(REAL_PATH "${CMAKE_CXX_COMPILER}" cxx_compiler_real_path)
      cmake_path(GET cxx_compiler_real_path PARENT_PATH cxx_compiler_dir)
      find_program(LLD_PATH NAMES ld.lld lld HINTS "${cxx_compiler_dir}")

      if(NOT LLD_PATH)
        # An lld older than clang may fail to link clang's LTO bitcode, so
        # name the package that matches the compiler (lld-22 for clang 22)
        string(REGEX MATCH "^[0-9]+" clang_major "${CMAKE_CXX_COMPILER_VERSION}")
        message(WARNING
          "lld not found in ${cxx_compiler_dir} or on the PATH; disabling IPO. "
          "Install the lld of the same version as clang, for example: "
          "sudo apt install lld-${clang_major}")
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

# ----------------------------- Main executable  -----------------------------
add_executable(cartogram src/main.cpp)
target_link_libraries(cartogram
  PRIVATE cartogram_lib
)
install(TARGETS cartogram RUNTIME DESTINATION bin)

# ----------------------------- Tests  -----------------------------
add_subdirectory(tests)
