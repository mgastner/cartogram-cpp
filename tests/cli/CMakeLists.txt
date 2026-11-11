# Tests of the command-line behavior of the cartogram binary. They carry the
# label "stress" because the stress-test job in CI runs the binary.

set(MAP_DIR "${CMAKE_SOURCE_DIR}/sample_data/usa_by_state_since_1959")

# With --make_csv, the program asks which property identifies the regions.
# If standard input ends before a valid answer, it must stop with
# EXIT_FAILURE and an error message instead of asking again forever.
add_test(
  NAME cli_make_csv_stops_at_end_of_input
  COMMAND sh -c "\"$0\" \"$1\" --make_csv < /dev/null 2> make_csv_stderr.txt; test $? -eq 1 && grep -q '^ERROR: ' make_csv_stderr.txt"
    "${CMAKE_BINARY_DIR}/cartogram"
    "${MAP_DIR}/usa_by_state_since_1959.geojson"
  WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
)
set_tests_properties(cli_make_csv_stops_at_end_of_input
  PROPERTIES
  LABELS "stress"
  TIMEOUT 60
)

# IDs with leading or trailing spaces in the GeoJSON must match the same IDs
# in the CSV, with or without the spaces (#310): no region may fall back to
# "no data", which the program reports as "ID ... in GeoJSON is not in CSV".
set(DATA_DIR "${CMAKE_CURRENT_SOURCE_DIR}/data")
foreach(csv_ids IN ITEMS padded trimmed)
  set(test_name cli_ids_with_spaces_match_${csv_ids}_csv_ids)
  set(work_dir "${CMAKE_CURRENT_BINARY_DIR}/${test_name}")
  file(MAKE_DIRECTORY "${work_dir}")
  add_test(
    NAME ${test_name}
    COMMAND sh -c "\"$0\" \"$1\" \"$2\" -P 100000000 > /dev/null 2> err.txt && ! grep -q 'not in CSV' err.txt"
      "${CMAKE_BINARY_DIR}/cartogram"
      "${DATA_DIR}/two_squares_padded_ids.geojson"
      "${DATA_DIR}/two_squares_${csv_ids}_ids.csv"
    WORKING_DIRECTORY "${work_dir}"
  )
  set_tests_properties(${test_name}
    PROPERTIES
    LABELS "stress"
    TIMEOUT 60
  )
endforeach()

# If the GeoJSON has both "Left" and "Left ", the CSV row "Left" belongs to
# the exact match "Left", whose target area must be the value in the CSV.
set(test_name cli_exact_id_match_takes_priority)
set(work_dir "${CMAKE_CURRENT_BINARY_DIR}/${test_name}")
file(MAKE_DIRECTORY "${work_dir}")
add_test(
  NAME ${test_name}
  COMMAND sh -c "\"$0\" \"$1\" \"$2\" -P 100000000 --export_preprocessed > /dev/null 2>&1 && grep -q '^Left,1.000000$' three_squares_exact_and_padded_ids_input_processed.csv"
    "${CMAKE_BINARY_DIR}/cartogram"
    "${DATA_DIR}/three_squares_exact_and_padded_ids.geojson"
    "${DATA_DIR}/three_squares_exact_and_padded_ids.csv"
  WORKING_DIRECTORY "${work_dir}"
)
set_tests_properties(${test_name}
  PROPERTIES
  LABELS "stress"
  TIMEOUT 60
)
