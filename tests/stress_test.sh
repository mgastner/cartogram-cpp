#!/usr/bin/env bash

# Start time, and other metadata
start_date=$(date '+%Y-%m-%d_%H-%M')

# Create results directory and change to it
results_dir="results_${start_date}"
mkdir -p "${results_dir}"
cd "${results_dir}"

results_file="log.txt"
successful_runs="successful_runs.txt"
failed_runs="failed_runs.txt"
tmp_file="tmp_file.txt"
SECONDS=0

# Arrays to store runtimes and map-csv pairs
runtimes=()
map_csv_pairs=()

# Add colors
red=1
green=2
yellow=3
blue=4
magenta=5
cyan=6
white=7
color() {
  tput setaf $1
  cat
  tput sgr0
}

# Flag to indicate if the script should skip the current test

# Function to handle SIGINT (Ctrl + C)
handle_sigint() {
  printf "\n\nSkipping current test...\n\n" | color $red
}

# Trap SIGINT (Ctrl + C)
trap handle_sigint SIGINT

# Parsing command line arguments
# From https://unix.stackexchange.com/questions/462515/
# check-command-line-arguments

verbose=0 # Initialize verbose flag to 0 (off)
flags=0 # Initialize a flag to test relevant cartogram-cpp flags to 0 (off)
cli_arr=("" "--plot_polygons" "--export_preprocessed" "--export_time_report" "--output_equal_area_map") # Set of flags to test if flags flag is enabled
cli=""    # Default CLI options
flags_arr=("No flags" "Plot polygons flag" "Export preprocessed flag" "Export time report flag" "Output equal area map flag")
total_arr=(0 0 0 0 0)
failed_arr=(0 0 0 0 0)

# Parsing command line arguments
if [ $# -eq 0 ]; then
  cli=""
else
  for arg in "$@"; do
    if [ "$arg" == "--verbose" ]; then
      verbose=1
      printf "\VERBOSE mode turned on.\n"
    elif [ "$arg" == "--flags" ]; then
      flags=1
      printf "\FLAGS mode turned on.\n"
    else
      cli="$cli $arg"
    fi
  done
fi

printf "\nWriting to ${results_file}\n"
printf "Tested on ${start_date}\n" >>"${results_file}"

if [ -n "$cli" ]; then
  printf "Passing" | tee -a "${results_file}"
  printf " ${cli} " | tee -a "${results_file}" | color $magenta
  printf "to cartogram\n\n" | tee -a "${results_file}"
else
  printf "No command line flags passed to /usr/bin/cartogram\n\n" | tee -a "${results_file}"
fi

# Turning on extended glob option
shopt -s extglob nocasematch

# Simple progress bar:
# Inspired from: https://github.com/pollev/bash_progress_bar

draw_progress_bar() {
  local terminal_width=$(tput cols)
  local progress_bar_width=0
  local max_progress_bar_width=50
  local padding=0
  if [[ "$terminal_width" -gt 60 ]]; then
    progress_bar_width=$(((terminal_width - 56) < max_progress_bar_width ? (terminal_width - 56) : max_progress_bar_width))
    padding=50
  else
    progress_bar_width=$((terminal_width - 6))
  fi
  local __percentage=$((10#$1))

  # Rescale the bar according to the progress bar width
  local __num_bar=$(($__percentage * $progress_bar_width / 100))

  # Add padding before the progress bar
  for p in $(seq 1 $padding); do printf " "; done

  # Draw progress bar
  for b in $(seq 1 $__num_bar); do printf "⬜"; done
  for s in $(seq 1 $(($progress_bar_width - $__num_bar))); do printf " "; done
  printf "$__percentage%% \r"
}

# Function to test map with csv
run_map() {

  printed=0
  integration_count=0
  max_area_err=""
  start=$SECONDS
  curr_cli=$cli
  # If world map then use -W flag
  if [[ "${country}" == world* ]]; then
    printf "World map detected. Passing -W flag.\n" | tee -a "${results_file_loc}" | color $magenta
    curr_cli="${cli} -W"
  fi

  printf "cartogram ${map} ${csv} ${curr_cli}\n\n"
  draw_progress_bar 0
  stdbuf -oL cartogram ${map} ${csv} ${curr_cli} 2>&1 |
    while read line; do
      # save to temp file
      echo $line >>${tmp_file}

      case "$line" in
      *"error"* | *"warning"* | *"invalid"*)
        printf "\n\n%s\n\n" "$line" | tee -a "${results_file_loc}" | color $red
        ;;
      *"progress: 0."*)
        draw_progress_bar ${line:12:2}
        ;;
      *"Integration number"*)
        integration_count=$((integration_count + 1))
        printf "%s\n" "$line" | color $red
        ;;
      *"New grid dimensions:"*)
        printf "%s\n" "$line" | color $blue
        ;;
      *"Max. area err:"*)
        max_area_err=$line
        printf "%s\n" "$line" | color $yellow
        ;;
      *"Area drift:"*)
        printf "%s\n" "$line" | color $magenta
        ;;
      *"progress: 1"*)
        if [[ "$printed" -eq 0 ]]; then
          draw_progress_bar 100
          printed=1
          printf "\n\n== Integration finished ==\nTotal integrations done: %d\n%s\n" "$integration_count" "$max_area_err" | tee -a "${results_file_loc}" | color $yellow
        fi
        ;;
      *" ms")
        printf "%s\n" "$line" | tee -a "${results_file_loc}" | color $cyan
        ;;
      # Check the output when using output_equal_area_map by line directly
      *"Writing"*)
        if [[ "${curr_cli}" == *output_equal_area_map* ]]; then
          printf "%s\n" "$line" | color $yellow
          draw_progress_bar 100
        fi
        ;;
      *)
        if [ $verbose -eq 1 ]; then
          printf "%s\n" "$line"
        fi
        ;;
      esac
    done
  end=$SECONDS
  runtime=$((end - start))
  printf "== Runtime ${runtime}s == \n" | tee -a "${results_file_loc}"

  # Checking for any errors, invalid geometry or unfinished integration
  # Since output_equal_area_map has no progress output and only has
  # a single line as output, check if line prints out with no errors
  if grep -qi "invalid" ${tmp_file} || \
     grep -qi "error" ${tmp_file} || \
     { ! grep -Fxq "Progress: 1" ${tmp_file} && [[ "${curr_cli}" != *output_equal_area_map* ]]; } || \
     { [[ "${curr_cli}" == *output_equal_area_map* ]] && ! grep -qi "writing" ${tmp_file}; }; then
    printf "== FAILED ==\n" | tee -a "${results_file_loc}" | color $red
    printf "cartogram ${map} ${csv} ${cli}" | tee -a ${failed_runs} | pbcopy
    printf "\n" >>${failed_runs}
    # Printing country to failed_tmp.txt
    if [[ "${failed}" -eq 0 ]] || ! grep -qi "${country}" failed_tmp.txt; then
      printf "\n${country}\n" >>failed_tmp.txt
    fi

    # Increasing failed counter
    failed_arr[i]=$((failed_arr[i] + 1))

    # Saving data to file if FAILED
    map_wo_ext=${map_file_name%.*json}
    csv_wo_ext=${csv_name%.csv}
    err_file="DNF-${map_wo_ext}-with-${csv_wo_ext}.txt"
    printf " - with ${csv_name} in ${runtime}s\n" >>failed_tmp.txt
    printf "cartogram %s %s %s\n" "$map" "$csv" "$cli" >>failed_tmp.txt
    printf "Full output saved to ${err_file}\n" | tee -a "${results_file_loc}"
    mv ${tmp_file} ${err_file}

    # If no errors, and integration finished, pass
  else
    printf "== PASSED ==\n" | tee -a "${results_file_loc}" | color $green
    printf "cartogram ${map} ${csv} ${cli}\n" >>${successful_runs}

    # Store the runtime for the current map
    runtimes+=($runtime)
    # Mention the flag used in the map csv pair
    # If no flag is used, print as "no flags"
    if [ -z "$cli" ]; then
      curr_cli="no flags"
    fi
    map_csv_pairs+=("${map_file_name} with ${csv_name} and ${curr_cli}")

  fi

  # Empty temporary file
  >${tmp_file}

  # Printing new line
  printf "\n" | tee -a "${results_file_loc}"
}

# Counting number of tests
countries=0

# If flags enabled, change the sample_data to its appropriate location
# Store results_file as one outside of subdirectories instead of splitting it
if [ $flags -eq 1 ]; then
  sample_data_folder="../../../sample_data/*"
  results_file_loc="../${results_file}"
# If flags not enabled, use single cli specified by user
else
  cli_arr=("$cli")
  passed_arr=(0)
  failed_arr=(0)
  sample_data_folder="../../sample_data/*"
  results_file_loc="${results_file}"
fi

# Loop through cli array
# Only loops once if flags not enabled
for i in "${!cli_arr[@]}"; do
  cli=${cli_arr[$i]}
  # If flags enabled, create and change to subdirectory for each flag
  if [ $flags -eq 1 ]; then
    if [ -n "$cli" ]; then
      dir_name="${cli:2}"
      mkdir -p "${dir_name}"
      cd "${dir_name}"
    else
      mkdir -p "no_flags"
      cd "no_flags"
    fi
  fi

  # Iterating through folders in sample_data folder
  for folder in ${sample_data_folder}; do
    if [[ -d "${folder}" && "${folder}" != *"sandbox"* ]]; then
      if [ $i -eq 0 ]; then
        countries=$((countries + 1))
      fi
      country=${folder##*/}
      printf " -------- Testing ${country}\n\n" | tee -a "${results_file_loc}" | color $magenta

      # Iterating through maps (GeoJSON(s) or JSON(s)) in country's folder
      for map in ${folder}/*.*json; do
        map_file_name=${map##*/}
        printf "Trying ${map_file_name}\n" | tee -a "${results_file_loc}" | color $blue

        # Iterating through visual variable files (CSV(s)) in country's folder
        for csv in ${folder}/*.csv; do
          total_arr[i]=$((total_arr[i] + 1))
          csv_name=${csv##*/}
          printf " - with ${csv_name}\n\n" | tee -a "${results_file_loc}" | color $cyan

          # Running cartogram, and timing it in seconds
          run_map

        done
      done
      printf " -------- Finished testing ${country}.\n\n\n\n\n\n" | tee -a "${results_file_loc}" | color $magenta
    fi
  done

  # Removing temporary files
  rm ${tmp_file}
  # Change directory to outside flag subdirectory once done with loop
  if [ $flags -eq 1 ]; then
    cd ..
  fi

done

# Summary report
printf "===== Finished testing all countries. =====\n\n" | tee -a "${results_file}" | color $magenta

duration="$(($SECONDS / 60))m $(($SECONDS % 60))s"
total_runtime=0
max_runtime=0
max_runtime_map=""
median_runtime=0

# Calculate total runtime, max runtime, and find the map with max runtime
for i in "${!runtimes[@]}"; do
  total_runtime=$((total_runtime + runtimes[i]))
  if ((runtimes[i] > max_runtime)); then
    max_runtime=${runtimes[i]}
    max_runtime_map=${map_csv_pairs[i]}
  fi
done

total_tests=0
for num in "${total_arr[@]}"; do
  total_tests=$((total_tests + num))
done

# Calculate average runtime
average_runtime=$((total_runtime / total_tests))

# Calculate median runtime
sorted_runtimes=($(printf '%s\n' "${runtimes[@]}" | sort -n))
if ((total_tests % 2 == 0)); then
  median_runtime=$(((sorted_runtimes[total_tests / 2 - 1] + sorted_runtimes[total_tests / 2]) / 2))
else
  median_runtime=${sorted_runtimes[total_tests / 2]}
fi

# Calculate the top quartile runtime maps
# top_quartile_count=$((total_tests / 4))
# if (( top_quartile_count == 0 )); then
#   top_quartile_count=1
# fi
# sorted_indices=($(printf '%s\n' "${!runtimes[@]}" | sort -rn -k1,1 -t ' '))
# top_quartile_maps=()
# for i in $(seq 0 $((top_quartile_count - 1))); do
#   idx=${sorted_indices[$i]}
#   top_quartile_maps+=("${runtimes[$idx]}s - ${map_csv_pairs[$idx]}")
# done

printf "Finished ${total_tests} tests on ${countries} countries in ${duration}.\n"
printf "Total runtime: ${total_runtime}s\n" | tee -a "${results_file}"
printf "Average runtime: ${average_runtime}s\n" | tee -a "${results_file}"
printf "Max runtime: ${max_runtime}s (Map: ${max_runtime_map})\n" | tee -a "${results_file}"
printf "Median runtime: ${median_runtime}s\n" | tee -a "${results_file}"

# printf "Top quartile maps with longest runtime:\n" | tee -a "${results_file}"
# for map in "${top_quartile_maps[@]}"; do
#   printf " - ${map}\n" | tee -a "${results_file}"
# done

total_failed=0
for num in "${failed_arr[@]}"; do
  total_failed=$((total_failed + num))
done

if [ $flags -eq 1 ]; then
  for i in "${!total_arr[@]}"; do
    failed=${failed_arr[$i]}
    total=${total_arr[$i]}

    printf "${flags_arr[$i]}:\n" | tee -a "${results_file}"
    failed_per=$((100 * $failed / $total))
    printf "Passed [$((total - failed))/${total}] | $((100 - failed_per))%% \n" | tee -a "${results_file}" | color $green
    printf "Failed [${failed}/${total}] | ${failed_per}%% \n" | tee -a "${results_file}" | color $red
  done
  printf "Combined tests:\n" | tee -a "${results_file}"
fi

failed_per=$((100 * $total_failed / $total_tests))
printf "Passed [$((total_tests - total_failed))/${total_tests}] | $((100 - failed_per))%% \n" | tee -a "${results_file}" | color $green
printf "Failed [${total_failed}/${total_tests}] | ${failed_per}%% \n" | tee -a "${results_file}" | color $red

# Checking if any tests failed
if [[ "${total_failed}" -gt 0 ]]; then

  # Showing failed tests and removing temporary file
  printf "\nFailed tests:\n" | tee -a "${results_file}" | color $magenta

  find . -type f -name "failed_tmp.txt" | while IFS= read -r file; do
    while IFS= read -r line; do
      if [[ $(echo $line | wc -w) -eq 1 ]]; then
        printf "$line\n" | tee -a "${results_file}" | color $red
      elif [[ $line == " - "* ]]; then
        printf "$line\n" | tee -a "${results_file}" | color $yellow
      else
        printf "$line\n" | tee -a "${results_file}"
      fi
    done < "$file"
    rm "$file"
  done
  printf "\n" | tee -a "${results_file}"
fi

printf "Results saved to ${results_dir}\n\n"

# Prompting for file deletion
printf "Clear ALL *.geojson, *.csv, *.svg, *.png and *.ps files in directory? [y/N]: " | color $yellow
read -t 10 to_clear
if [[ $? -ne 0 ]]; then
  printf "\nNo response received! " | color $red
elif [ "$to_clear" == "y" ]; then
  # Clear the files in subdirectories as well if flags are enabled
  find . -name "*.geojson" -type f -delete
  find . -name "*.csv" -type f -delete
  find . -name "*.svg" -type f -delete
  find . -name "*.png" -type f -delete
  find . -name "*.ps" -type f -delete
  printf "All *.geojson, *.csv, *.svg, *.png and *.ps files deleted.\n" | color $red
  exit ${total_failed}
fi

printf "Files not cleared.\n" | color $green
exit ${total_failed}
