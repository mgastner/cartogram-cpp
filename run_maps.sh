#!/bin/bash

# Save all the flags and arguments passed to the script
user_args=("$@")

# Change to the directory containing the sample_data folder and exit if it fails
pushd "$(dirname "$0")" >/dev/null || exit
cd sample_data || exit

# Initialize success and fail counters
success_count=0
fail_count=0

print_separator() {
  echo "----------------------------------------------------------------------------------------------------"
}

# Iterate through the directories inside sample_data
for map_dir in */; do
  # Check if there's a .geojson file in the directory
  geojson_file=$(find "$map_dir" -type f -iname "*.geojson" | head -n 1)
  if [ -z "$geojson_file" ]; then
    continue
  fi

  # Find all .csv files in the directory
  csv_files=$(find "$map_dir" -type f -iname "*.csv")
  if [ -z "$csv_files" ]; then
    continue
  fi

  # Process each CSV file in the directory iteratively
  for csv_file in $csv_files; do
    print_separator

    # Remove trailing slash from map_dir
    map_name=$(basename "$map_dir")
    # Only keep the part after the last / in the CSV file path
    data_file=$(basename "$csv_file")

    echo "Map: $map_name"
    echo "Data: $data_file"

    cmd=(../build/Release/cartogram
      "../sample_data/$geojson_file"
      "../sample_data/$csv_file")

    [[ $map_name == *world* ]] && cmd+=(--world)
    cmd+=("${user_args[@]}")

    printf 'Command: '
    printf '%q ' "${cmd[@]}"
    echo
    echo "Current folder: $(pwd)/$map_dir"

    # Run the command and capture stderr output while timing it
    SECONDS=0
    stderr_output=$("${cmd[@]}" 2>&1 >/dev/null)

    # Extract the last occurrences of "Integration number" and "Total time:" from stderr_output
    integration_line=$(echo "$stderr_output" | grep "Integration number" | tail -n 1)
    total_time_line=$(echo "$stderr_output" | grep "Total time:" | tail -n 1)

    # Check for errors first; then determine success based on "Total time"
    if echo "$stderr_output" | grep -q "ERROR"; then
      fail_count=$((fail_count + 1))
      echo "Result: Fail"
      echo ""
      echo "Integration: $integration_line"
    elif echo "$stderr_output" | grep -q "Total time"; then
      success_count=$((success_count + 1))
      echo "Result: Success"
      echo ""
      echo "$integration_line"
      echo "$total_time_line"
    else
      fail_count=$((fail_count + 1))
      echo "Result: Fail"
      echo ""
      echo "Integration: $integration_line"
    fi

    print_separator
    echo ""
  done
done

# Print the final success and fail counts
echo "Final Results:"
echo "Success count: $success_count"
echo "Fail count: $fail_count"

# Return to the original directory and exit if it fails
popd >/dev/null || exit
