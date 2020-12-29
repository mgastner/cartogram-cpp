if (MSVC)
    # warning level 4
    add_compile_options(/W4)
else()
    # lots of warnings
    add_compile_options(-Wall -Wextra -pedantic)
endif()

# code for findFFTW
configure_file(downloadFindFFTW.cmake.in findFFTW-download/CMakeLists.txt)
execute_process(COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" .
        RESULT_VARIABLE result
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/findFFTW-download )
if(result)
    message(FATAL_ERROR "CMake step for findFFTW failed: ${result}")
    else()
    message("CMake step for findFFTW completed (${result}).")
endif()
execute_process(COMMAND ${CMAKE_COMMAND} --build .
        RESULT_VARIABLE result
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/findFFTW-download )
if(result)
    message(FATAL_ERROR "Build step for findFFTW failed: ${result}")
endif()

set(findFFTW_DIR ${CMAKE_CURRENT_BINARY_DIR}/findFFTW-src)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${findFFTW_DIR}")

# As per https://github.com/CGAL/cgal/wiki/How-to-use-CGAL-with-CMake-or-your-own-build-system
cmake_minimum_required(VERSION 3.1)
project(cartogram LANGUAGES CXX)

# set(PLATFORM_SPECIFIC_LIBS "-lpthread")
set(CMAKE_CXX_STANDARD 20)
# add_compile_definitions(_GLIBCXX_USE_CXX11_ABI=0)
set(THREADS_PREFER_PTHREAD_FLAG ON)
set(CMAKE_BUILD_TYPE Release)
set(CMAKE_CXX_COMPILER "g++-10")
set(CMAKE_C_COMPILER "gcc-10")

find_package(CGAL)
# include system to boost?
find_package(Boost 1.70 COMPONENTS program_options REQUIRED)
find_package(FFTW)
find_package(Threads REQUIRED)
find_package(OpenMP REQUIRED)

# include_directories(${Boost_INCLUDE_DIR})
# link_directories(${Boost_LIBRARY_DIR})

add_executable(cartogram ../cartogram_generator/main.cpp
              ../cartogram_generator/geo_div.cpp
          		../cartogram_generator/map_state.cpp
          		../cartogram_generator/read_csv.cpp
          		../cartogram_generator/rescale_map.cpp
          		../cartogram_generator/write_eps.cpp
              # ../cartogram_generator/topological_checks.cpp
              # ../cartogram_generator/helper_functions.cpp
          		../cartogram_generator/read_geojson.cpp
          		../cartogram_generator/fill_with_density.cpp
              ../cartogram_generator/blur_density.cpp
              ../cartogram_generator/flatten_density.cpp
            ../cartogram_generator/ft_real_2d.cpp)

# Boost::program_options ${Boost_LIBRARIES}
target_link_libraries(cartogram
                      ${Boost_LIBRARIES}
                      # ${Boost_SYSTEM_LIBRARY}
                      # Boost::program_options
                      OpenMP::OpenMP_CXX
                      CGAL::CGAL
                      ${FFTW_LIBRARIES}
                      ${CMAKE_THREAD_LIBS_INIT})
