find_package(Python3 COMPONENTS Interpreter REQUIRED)

set(MAP_ROOT "${CMAKE_SOURCE_DIR}/sample_data")
file(GLOB MAP_DIRS
  RELATIVE "${MAP_ROOT}"
  CONFIGURE_DEPENDS
  "${MAP_ROOT}/*"
)

foreach(map_dir IN LISTS MAP_DIRS)
  if(IS_DIRECTORY "${MAP_ROOT}/${map_dir}")
    add_test(NAME fuzz_${map_dir}
      COMMAND ${Python3_EXECUTABLE}
      ${CMAKE_SOURCE_DIR}/tests/fuzzer/fuzzer.py
      ${MAP_ROOT}
      ${CMAKE_BINARY_DIR}/cartogram
      --max-rows 100
      --runs 5
      --out-dir ${CMAKE_BINARY_DIR}/fuzzer_csvs
      --map ${map_dir}
    )
    set_tests_properties(fuzz_${map_dir} PROPERTIES
      LABELS "fuzzer"
      TIMEOUT 120
      SKIP_RETURN_CODE 77
    )
  endif()
endforeach()
