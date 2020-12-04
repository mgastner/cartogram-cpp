
if(NOT "/mnt/d/projects/github/cartogram_cpp_master/cartogram_generator/findFFTW-download/findFFTW_download-prefix/src/findFFTW_download-stamp/findFFTW_download-gitinfo.txt" IS_NEWER_THAN "/mnt/d/projects/github/cartogram_cpp_master/cartogram_generator/findFFTW-download/findFFTW_download-prefix/src/findFFTW_download-stamp/findFFTW_download-gitclone-lastrun.txt")
  message(STATUS "Avoiding repeated git clone, stamp file is up to date: '/mnt/d/projects/github/cartogram_cpp_master/cartogram_generator/findFFTW-download/findFFTW_download-prefix/src/findFFTW_download-stamp/findFFTW_download-gitclone-lastrun.txt'")
  return()
endif()

execute_process(
  COMMAND ${CMAKE_COMMAND} -E remove_directory "/mnt/d/projects/github/cartogram_cpp_master/cartogram_generator/findFFTW-src"
  RESULT_VARIABLE error_code
  )
if(error_code)
  message(FATAL_ERROR "Failed to remove directory: '/mnt/d/projects/github/cartogram_cpp_master/cartogram_generator/findFFTW-src'")
endif()

# try the clone 3 times in case there is an odd git clone issue
set(error_code 1)
set(number_of_tries 0)
while(error_code AND number_of_tries LESS 3)
  execute_process(
    COMMAND "/usr/bin/git"  clone --no-checkout "https://github.com/egpbos/findfftw.git" "findFFTW-src"
    WORKING_DIRECTORY "/mnt/d/projects/github/cartogram_cpp_master/cartogram_generator"
    RESULT_VARIABLE error_code
    )
  math(EXPR number_of_tries "${number_of_tries} + 1")
endwhile()
if(number_of_tries GREATER 1)
  message(STATUS "Had to git clone more than once:
          ${number_of_tries} times.")
endif()
if(error_code)
  message(FATAL_ERROR "Failed to clone repository: 'https://github.com/egpbos/findfftw.git'")
endif()

execute_process(
  COMMAND "/usr/bin/git"  checkout master --
  WORKING_DIRECTORY "/mnt/d/projects/github/cartogram_cpp_master/cartogram_generator/findFFTW-src"
  RESULT_VARIABLE error_code
  )
if(error_code)
  message(FATAL_ERROR "Failed to checkout tag: 'master'")
endif()

set(init_submodules TRUE)
if(init_submodules)
  execute_process(
    COMMAND "/usr/bin/git"  submodule update --recursive --init 
    WORKING_DIRECTORY "/mnt/d/projects/github/cartogram_cpp_master/cartogram_generator/findFFTW-src"
    RESULT_VARIABLE error_code
    )
endif()
if(error_code)
  message(FATAL_ERROR "Failed to update submodules in: '/mnt/d/projects/github/cartogram_cpp_master/cartogram_generator/findFFTW-src'")
endif()

# Complete success, update the script-last-run stamp file:
#
execute_process(
  COMMAND ${CMAKE_COMMAND} -E copy
    "/mnt/d/projects/github/cartogram_cpp_master/cartogram_generator/findFFTW-download/findFFTW_download-prefix/src/findFFTW_download-stamp/findFFTW_download-gitinfo.txt"
    "/mnt/d/projects/github/cartogram_cpp_master/cartogram_generator/findFFTW-download/findFFTW_download-prefix/src/findFFTW_download-stamp/findFFTW_download-gitclone-lastrun.txt"
  RESULT_VARIABLE error_code
  )
if(error_code)
  message(FATAL_ERROR "Failed to copy script-last-run stamp file: '/mnt/d/projects/github/cartogram_cpp_master/cartogram_generator/findFFTW-download/findFFTW_download-prefix/src/findFFTW_download-stamp/findFFTW_download-gitclone-lastrun.txt'")
endif()

