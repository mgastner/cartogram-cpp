#!/bin/bash

# Change to the directory containing the sample_data folder
pushd "$(dirname "$0")" > /dev/null
cd sample_data

# Initialize success and fail counters
success_count=0
fail_count=0

# Iterate through the directories inside sample_data
for map_dir in */; do
  # Check if there's a .geojson file in the directory
  geojson_file=$(find "$map_dir" -type f -iname "*.geojson" | head -n 1)
  
  if [ -z "$geojson_file" ]; then
    continue
  fi

  # Find the .csv file in the directory
  csv_file=$(find "$map_dir" -type f -iname "*.csv" | head -n 1)

  # Print the current map being processed
  echo "Processing map: $map_dir"

  # Construct the command to run
  command="./bin/cartogram \"${geojson_file}\" \"${csv_file}\" -Qrts"

  # Print the exact command being run
  # echo "Running command: $command"
  
  # add ./sample_data/ to the command before geojson_file and csv_file
  command_2="./bin/cartogram \"./sample_data/${geojson_file}\" \"./sample_data/${csv_file}\" -Qrts"
  
  # Print the exact command being run
  echo "Running command: $command_2"
  
  # Run the command and store stderr output in a variable, while timing the command
  SECONDS=0
  stderr_output=$(eval "../$command" 2>&1 >/dev/null)
  elapsed_time=$(awk "BEGIN {print $SECONDS * 1000}")

  # Check if "Total time" is in the stderr output
  if echo "$stderr_output" | grep -q "Total Time"; then
    success_count=$((success_count + 1))
    echo "Result: Success"
  else
    fail_count=$((fail_count + 1))
    echo "Result: Fail"
  fi

  # Print the time it took to process the map
  echo "Time taken: ${elapsed_time}ms"
  echo ""
done

# Print the final success and fail counts
echo "Success count: $success_count"
echo "Fail count: $fail_count"

# Return to the original directory
popd > /dev/null