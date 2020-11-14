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
set(CMAKE_BUILD_TYPE Release)

# As per https://github.com/CGAL/cgal/wiki/How-to-use-CGAL-with-CMake-or-your-own-build-system

# set(FFTW_INCLUDE_DIRS "/usr")

cmake_minimum_required(VERSION 3.1)
project(cartogram)

# set(PLATFORM_SPECIFIC_LIBS "-lpthread")
set(CMAKE_CXX_FLAGS "-std=c++20 -pthread")

find_package(CGAL)
find_package(Boost 1.40 COMPONENTS program_options REQUIRED)
find_package(FFTW)

add_executable(cartogram main.cpp
    geo_div.cpp
		map_state.cpp
		read_csv.cpp
		rescale_map.cpp
		write_eps.cpp
		read_geojson.cpp
		fill_with_density.cpp
    blur_density.cpp
    flatten_density.cpp
		ft_real_2d.cpp)


target_link_libraries(cartogram CGAL::CGAL Boost::program_options fftw3)
