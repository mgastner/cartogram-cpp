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
