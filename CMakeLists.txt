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
    set(lto_compile_options
      -flto=full -fwhole-program-vtables -fvirtual-function-elimination)

    if(ipo_ok)
      # Find the ld.lld that clang would run, make clang run exactly that
      # one, and check with a test link that it links clang's bitcode. The
      # test program is compiled for the same CPU as the library.
      include("${CMAKE_CURRENT_LIST_DIR}/cmake/LtoLinker.cmake")
      get_target_property(cpu_options cartogram_lib COMPILE_OPTIONS)
      list(FILTER cpu_options INCLUDE REGEX "^-m(arch|tune)=")
      cartogram_find_lto_linker(lto_link_options lto_linker
        ${cpu_options} ${lto_compile_options})

      if(NOT lto_link_options)
        set(ipo_ok FALSE)
      endif()
    endif()
  endif()

  if(ipo_ok)
    set_property(TARGET cartogram_lib PROPERTY INTERPROCEDURAL_OPTIMIZATION ON)
    target_link_options(cartogram_lib PUBLIC -flto=full)

    if(CMAKE_CXX_COMPILER_ID STREQUAL "Clang" AND CMAKE_SYSTEM_NAME STREQUAL "Linux")
      message(STATUS "Using LLD with full LTO + whole-program vtables: ${lto_linker}")
      target_link_options(cartogram_lib PUBLIC ${lto_link_options})
      target_compile_options(cartogram_lib PUBLIC ${lto_compile_options})
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
