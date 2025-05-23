set(MAP_ROOT "${CMAKE_SOURCE_DIR}/sample_data")

file(GLOB CHILDREN
  RELATIVE "${MAP_ROOT}"
  CONFIGURE_DEPENDS
  "${MAP_ROOT}/*"
)

foreach(map_dir IN LISTS CHILDREN)
  if(NOT IS_DIRECTORY "${MAP_ROOT}/${map_dir}")
    continue()
  endif()

  file(GLOB GEO
    RELATIVE "${MAP_ROOT}/${map_dir}"
    "${MAP_ROOT}/${map_dir}/*.geojson")

  if(GEO)
    list(GET GEO 0 geofile)
  else()
    message(WARNING "tests/stress: no .geojson in ${map_dir}, skipping")
    continue()
  endif()

  file(GLOB CSVs
    RELATIVE "${MAP_ROOT}/${map_dir}"
    "${MAP_ROOT}/${map_dir}/*.csv")

  if(NOT CSVs)
    message(WARNING "tests/stress: no .csv in ${map_dir}, skipping")
    continue()
  endif()

  foreach(csv IN LISTS CSVs)
    string(REPLACE "." "_" test_name "stress_${map_dir}_${csv}")

    set(TEST_CMD
      "${CMAKE_BINARY_DIR}/cartogram"
      "${MAP_ROOT}/${map_dir}/${geofile}"
      "${MAP_ROOT}/${map_dir}/${csv}"
    )

    if(map_dir MATCHES "world")
      list(APPEND TEST_CMD "--world")
    endif()

    add_test(
      NAME ${test_name}
      COMMAND ${TEST_CMD}
    )
    set_tests_properties(${test_name} PROPERTIES LABELS "stress" TIMEOUT 120)
  endforeach()
endforeach()
