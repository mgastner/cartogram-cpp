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

  # Skip sandbox directory
  if(map_dir STREQUAL "sandbox")
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
    set(TEST_TIMEOUT 600)

    set_tests_properties(${test_name}
      PROPERTIES
      LABELS "stress"
      TIMEOUT ${TEST_TIMEOUT}
    )
  endforeach()
endforeach()

# Without --world, the world map by country used to crash with a segmentation
# fault in the simplification, while its triangulation used plain
# floating-point predicates. The test runs in its own folder, because its
# output has the same name as that of the test with --world.
set(world_dir "${MAP_ROOT}/world_by_country_since_2022")
set(without_world_dir "${CMAKE_CURRENT_BINARY_DIR}/without_world_option")
file(MAKE_DIRECTORY "${without_world_dir}")
add_test(
  NAME stress_world_by_country_since_2022_without_world_option
  COMMAND
    "${CMAKE_BINARY_DIR}/cartogram"
    "${world_dir}/world_by_country_since_2022.geojson"
    "${world_dir}/world_population_by_country_2010.csv"
  WORKING_DIRECTORY "${without_world_dir}"
)
set_tests_properties(stress_world_by_country_since_2022_without_world_option
  PROPERTIES
  LABELS "stress"
  TIMEOUT 600
)
