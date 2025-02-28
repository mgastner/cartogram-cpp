#!/usr/bin/env bash

# Initialize colors for output
red=1
green=2
yellow=3
color() {
  tput setaf "$1"
  cat
  tput sgr0
}

# Create or overwrite the output CSV file
output_csv="output.csv"
echo "Map Name,Min Density,Max Density,Difference" >"$output_csv"

# Function to run the cartogram command and capture densities
run_map() {
  map="$1"
  csv="$2"

  # Get the map name (filename without path)
  map_name=$(basename "$csv")

  # Display the command being run
  printf "Running command: cartogram %s %s\n\n" "$map" "$csv" | color "$yellow"

  # Initialize variables
  local min_density=""
  local max_density=""


  while IFS= read -r line; do
    echo "$line"
    case "$line" in
      *"Density min:"*)
        # Extract the value after "Density min:"
        min_density=${line#*Density min: }
        echo "Min Density: $min_density" | color "$green"
        ;;
      *"Density max:"*)
        # Extract the value after "Density max:"
        max_density=${line#*Density max: }
        echo "Max Density: $max_density" | color "$green"
        ;;
        *"Density mean:"*)
        # Extract the value after "Density mean:"
        mean_density=${line#*Density mean: }
        echo "mean Density: $mean_density" | color "$green"
        ;;
    esac
  done < <("../bin/cartogram" "${map}" "${csv}" 2>&1)
  
  ratio=$(echo "$max_density / $min_density" | bc)
  # Write to the CSV file even if values are not numeric
  echo "$map_name,$min_density,$max_density,$ratio" >>"$output_csv"
}

# Iterate through folders and files
for folder in ../sample_data/*; do
  if [[ -d "${folder}" && "${folder}" != *"sandbox"* ]]; then
    for map in "${folder}"/*.*json; do
      for csv in "${folder}"/*.csv; do
        run_map "${map}" "${csv}"
      done
    done
  fi
done

# Display the final CSV output
printf "\nOutput saved to %s\n" "$output_csv" | color "$green"
