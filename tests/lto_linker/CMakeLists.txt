# Small LTO build with the linker selection of the top-level CMakeLists.txt
# (cmake/LtoLinker.cmake), for check_lto_linker.sh. It is not part of the
# main build.
cmake_minimum_required(VERSION 3.30)
project(lto_linker_check LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

include(CheckIPOSupported)
check_ipo_supported(RESULT ipo_ok OUTPUT ipo_msg)
if(NOT ipo_ok)
  message(FATAL_ERROR "IPO/LTO not available: ${ipo_msg}")
endif()

include("${CMAKE_CURRENT_LIST_DIR}/../../cmake/LtoLinker.cmake")

add_library(shapes STATIC shapes.cpp)
target_include_directories(shapes PUBLIC "${CMAKE_CURRENT_LIST_DIR}")
add_executable(lto_linker_check main.cpp)
target_link_libraries(lto_linker_check PRIVATE shapes)

set(lto_compile_options
  -flto=full -fwhole-program-vtables -fvirtual-function-elimination)
cartogram_find_lto_linker(lto_link_options lto_linker ${lto_compile_options})

if(lto_link_options)
  set_property(TARGET shapes PROPERTY INTERPROCEDURAL_OPTIMIZATION ON)
  target_compile_options(shapes PUBLIC ${lto_compile_options})
  target_link_options(shapes PUBLIC -flto=full ${lto_link_options})
endif()

message(STATUS "LTO linker: ${lto_linker}")
