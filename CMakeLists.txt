cmake_minimum_required(VERSION 3.6)

project(cartogram LANGUAGES CXX)

# ========== Project Setup ==========
set(CMAKE_CXX_STANDARD 20)
set(THREADS_PREFER_PTHREAD_FLAG ON)
set(CMAKE_COLOR_DIAGNOSTICS ON)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# Default build type
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING "" FORCE)
endif()

# ========== Dependencies Setup ==========

# Boost
find_package(Boost REQUIRED COMPONENTS unit_test_framework)

# Matplot++
find_package(Matplot++)

# PkgConfig, fftw, and cairo
find_package(PkgConfig REQUIRED)
pkg_search_module(fftw REQUIRED fftw3 IMPORTED_TARGET)
pkg_search_module(cairo REQUIRED cairo IMPORTED_TARGET)

# ========== Source Files ==========
file(GLOB_RECURSE CARTOGRAM_SOURCES "src/*.cpp")
add_executable(cartogram ${CARTOGRAM_SOURCES})

# ========== Include Directories ==========
target_include_directories(cartogram PUBLIC
  ${PROJECT_SOURCE_DIR}/include
  ${Boost_INCLUDE_DIRS}
  PkgConfig::fftw
  PkgConfig::cairo
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
  PkgConfig::cairo
  $<$<BOOL:${Matplot++_FOUND}>:Matplot++::matplot>
)

# ========== Installation ==========
install(TARGETS cartogram DESTINATION bin)

# Enable CTest testing
enable_testing()

# Add test executable
file(GLOB_RECURSE TEST_FILES "tests/*.cpp")

# Include the source files from the src directory that are needed for testing
set(CARTOGRAM_TEST_SOURCES_FROM_SRC
  "src/misc/string_to_decimal_converter.cpp"

  # Add additional test sources from src here if necessary
)

# For each test file, create an executable and a test
foreach(TEST_FILE ${TEST_FILES})
  # Extract the filename without an extension to use as a test name
  get_filename_component(TEST_NAME ${TEST_FILE} NAME_WE)

  # Create an executable for each test file
  add_executable(${TEST_NAME} ${TEST_FILE} ${CARTOGRAM_TEST_SOURCES_FROM_SRC})

  # Include directories for the test executable
  target_include_directories(${TEST_NAME} PUBLIC
    ${PROJECT_SOURCE_DIR}/include
    ${Boost_INCLUDE_DIRS}
    PkgConfig::fftw
  )

  if(CMAKE_COMPILER_IS_GNUCXX OR CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    target_compile_options(${TEST_NAME} PRIVATE -isystem ${Boost_INCLUDE_DIRS})
    target_compile_options(${TEST_NAME} PRIVATE -ffp-contract=off)
  elseif(MSVC)
    target_compile_options(${TEST_NAME} PRIVATE /external:I ${Boost_INCLUDE_DIRS})
  endif()

  # Compiler warnings for the test executable
  target_compile_options(${TEST_NAME} PRIVATE -Wall -Wextra -pedantic -Wno-deprecated-declarations)

  # Linking Libraries for Test Executable
  target_link_libraries(${TEST_NAME}
    ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY}
    PkgConfig::fftw
  )

  # Register the executable as a test
  add_test(NAME ${TEST_NAME} COMMAND ${TEST_NAME})
endforeach()

# Uninstall target
add_custom_target("uninstall")
add_custom_command(
  TARGET "uninstall"
  POST_BUILD
  COMMENT "Uninstalling cartogram..."
  COMMAND xargs rm -vf < install_manifest.txt || echo "Nothing in install_manifest.txt to be uninstalled!"
)