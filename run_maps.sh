#!/bin/bash

# Change to the directory containing the sample_data folder and exit if it fails
pushd "$(dirname "$0")" >/dev/null || exit
cd sample_data || exit

# Initialize success and fail counters
success_count=0
fail_count=0

# Sum of all command runtimes (in seconds)
total_secs=0

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

    extra=""
    if [[ $map_name == *world* ]]; then
      extra="--world"
    fi

    command="../build/Release/cartogram \"../sample_data/${geojson_file}\" \"../sample_data/${csv_file}\" $extra"
    echo "Command: $command"

    # Show the full path from the script location for the current map folder
    current_folder="$(pwd)/$map_dir"
    echo "Current folder: $current_folder"

    # Run the command and capture stderr output while timing it
    SECONDS=0
    stderr_output=$(eval "$command" 2>&1 >/dev/null)
    elapsed_secs=$SECONDS
    total_secs=$((total_secs + elapsed_secs))

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

# Pretty-print total_secs as MM:SS
printf -v total_ms "%02d:%02d" $(((total_secs%3600)/60)) $((total_secs%60))
echo "Total time: ${total_ms} (${total_secs}s)"

# Return to the original directory and exit if it fails
popd >/dev/null || exit
