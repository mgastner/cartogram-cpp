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
# in the CSV, with or without the spaces, and with or without quotes (#310).
# The CSV reader removes the spaces even inside quotes, so each CSV ID
# matches the only GeoJSON ID that is equal to it after trimming. No region
# may fall back to "no data", which the program reports as "ID ... in
# GeoJSON is not in CSV", and each region must get its target area from the
# CSV.
set(DATA_DIR "${CMAKE_CURRENT_SOURCE_DIR}/data")
foreach(csv_ids IN ITEMS padded quoted trimmed)
  set(test_name cli_ids_with_spaces_match_${csv_ids}_csv_ids)
  set(work_dir "${CMAKE_CURRENT_BINARY_DIR}/${test_name}")
  file(MAKE_DIRECTORY "${work_dir}")
  add_test(
    NAME ${test_name}
    COMMAND sh -c "\"$0\" \"$1\" \"$2\" -P 100000000 --export_preprocessed > /dev/null 2> err.txt && ! grep -q 'not in CSV' err.txt && grep -Fqx 'Left  ,1.000000' \"$3\" && grep -Fqx ' Right,2.000000' \"$3\""
      "${CMAKE_BINARY_DIR}/cartogram"
      "${DATA_DIR}/two_squares_padded_ids.geojson"
      "${DATA_DIR}/two_squares_${csv_ids}_ids.csv"
      "two_squares_${csv_ids}_ids_input_processed.csv"
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

# Runs the program on a GeoJSON file and a CSV file in data/. The test passes
# if the program stops with EXIT_FAILURE and its standard error contains
# expected_text.
function(add_id_error_test test_name geojson csv expected_text)
  set(work_dir "${CMAKE_CURRENT_BINARY_DIR}/${test_name}")
  file(MAKE_DIRECTORY "${work_dir}")
  add_test(
    NAME ${test_name}
    COMMAND sh -c "\"$0\" \"$1\" \"$2\" -P 100000000 > /dev/null 2> err.txt; test $? -eq 1 && grep -Fq -- \"$3\" err.txt"
      "${CMAKE_BINARY_DIR}/cartogram"
      "${DATA_DIR}/${geojson}"
      "${DATA_DIR}/${csv}"
      "${expected_text}"
    WORKING_DIRECTORY "${work_dir}"
  )
  set_tests_properties(${test_name}
    PROPERTIES
    LABELS "stress"
    TIMEOUT 60
  )
endfunction()

# If the GeoJSON has "Left " and " Left" but not "Left", the CSV row "Left"
# could belong to either region. The program must stop instead of choosing
# by the order of the features, with the same message for both orders.
set(ambiguity_message "ERROR: ID \"Left\" in CSV matches 2 IDs in GeoJSON if leading and trailing whitespace is ignored: \" Left\", \"Left \".")
add_id_error_test(
  cli_ambiguous_id_match_fails
  three_squares_ids_equal_after_trimming.geojson
  three_squares_ids_equal_after_trimming.csv
  "${ambiguity_message}"
)
add_id_error_test(
  cli_ambiguous_id_match_fails_for_reversed_features
  three_squares_ids_equal_after_trimming_reversed.geojson
  three_squares_ids_equal_after_trimming.csv
  "${ambiguity_message}"
)

# Two CSV rows must not both give data to one region, whether their IDs are
# equal ("Left" twice) or equal after trimming ("Left" and "Left  ").
add_id_error_test(
  cli_duplicate_csv_ids_fail
  three_squares_exact_and_padded_ids.geojson
  three_squares_duplicate_ids.csv
  "ERROR: More than one row in CSV matches ID \"Left\" in GeoJSON"
)
add_id_error_test(
  cli_csv_ids_equal_after_trimming_fail
  two_squares_padded_ids.geojson
  two_squares_duplicate_ids.csv
  "ERROR: More than one row in CSV matches ID \"Left  \" in GeoJSON"
)

# A CSV column with a repeated ID ("Left" twice in the column "name") can
# match a GeoJSON property whose IDs are equal after trimming ("Left " and
# " Left"). The program must still choose the column "code", in which every
# ID appears once, and give the target areas to the regions by their codes.
set(test_name cli_id_column_without_repeated_ids_takes_priority)
set(work_dir "${CMAKE_CURRENT_BINARY_DIR}/${test_name}")
file(MAKE_DIRECTORY "${work_dir}")
add_test(
  NAME ${test_name}
  COMMAND sh -c "\"$0\" \"$1\" \"$2\" -P 100000000 --export_preprocessed > /dev/null 2>&1 && grep -Fqx 'L1,1.000000' three_squares_with_codes_input_processed.csv && grep -Fqx 'L2,3.000000' three_squares_with_codes_input_processed.csv"
    "${CMAKE_BINARY_DIR}/cartogram"
    "${DATA_DIR}/three_squares_with_codes.geojson"
    "${DATA_DIR}/three_squares_with_codes.csv"
  WORKING_DIRECTORY "${work_dir}"
)
set_tests_properties(${test_name}
  PROPERTIES
  LABELS "stress"
  TIMEOUT 60
)
