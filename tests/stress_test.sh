#!/usr/bin/env bash

# Start time, and other metadata
start_date=$(date '+%Y-%m-%d_%H-%M-%S')
results_file="results_${start_date}.txt"
tmp_file="tmp_file.txt"
SECONDS=0
printf "\nWriting to ${results_file}\n"
printf "Tested on ${start_date}\n" >> "${results_file}"

# Add colors
red=1; green=2; yellow=3; blue=4; magenta=5; cyan=6; white=7
color() { tput setaf $1; cat; tput sgr0; }

# Parsing command line arguments
# From https://unix.stackexchange.com/questions/462515/
# check-command-line-arguments
if [ $# -eq 0 ] || [ -z "$1" ]; then
cli="-est"
else
  cli="$1"
fi
printf "Passing" | tee -a "${results_file}"
printf " ${cli} " | tee -a "${results_file}" | color $magenta
printf "to cartogram\n\n" | tee -a "${results_file}"

# Turning on extended glob option
shopt -s extglob nocasematch

# Simple loading spinner, adapted from:
# https://stackoverflow.com/questions/238073/how-to-add-a-progress-bar-to-a-shell-script/4276265#4276265

sp='/-\|'
spin()
{
  for run in {1..3}; do
    printf '\b%.1s' "$sp"
    sleep 0.33
    sp=${sp#?}${sp%???}
  done
}

# Function to test map with csv
run_map()
{
  start=$SECONDS
  curr_runtime=$SECONDS
  cartogram ${map} ${csv} ${cli} 2>&1 |
    while read line
  do
    # save to temp file
    echo $line >> ${tmp_file}

    # display warnings, errors etc.
    if [[ $line =~ "error"  || $line =~ "warning" || $line =~ "invalid" ]]; then
    printf "\b$line\n\n" | tee -a "${results_file}" | color $red
    fi

    if [[ $line =~ "progress: 1" ]]; then
      printf "\b== Integration finished ==\n" | tee -a "${results_file}"
    fi

    # check if integration finished
    if [[ $line =~ "progress" ]]; then
    printf "\b - $line, in $((SECONDS-start))s\n" | tee -a "${results_file}"
    fi

    # code for spinner
    if [[ $(( SECONDS - $curr_runtime)) -eq 1 ]]; then
    spin &
    curr_runtime=$SECONDS
    fi

    done
  end=${SECONDS}
  runtime=$((end-start))
  printf "\b == Runtime ${runtime}s == \n" | tee -a "${results_file}"

  # Checking for any errors, invalid geometry or unfinished integration
  if grep -qi "invalid" ${tmp_file} || grep -qi "error" ${tmp_file} || ! grep -Fxq "Progress: 1" ${tmp_file} ; then
  printf "\b== FAILED ==\n" | tee -a "${results_file}" | color $red

  # Printing country to failed_tmp.txt
  if [[ "${failed}" -eq 0 ]] || ! grep -qi "${country}" failed_tmp.txt; then
  printf "\n${country}\n" >> failed_tmp.txt
  fi

  # Increasing failed counter
  failed=$((failed+1))

  # Saving data to file if FAILED
  map_wo_ext=${map_file_name%.*json}
  csv_wo_ext=${csv_name%.csv}
  err_file="results_${start_date}-${map_wo_ext}-${csv_wo_ext}.txt"
  printf " - ${map_file_name} with ${csv_name}\n" >> failed_tmp.txt
  printf "Full output saved to ${err_file}\n" | tee -a "${results_file}"
  mv ${tmp_file} ${err_file}

  # If no errors, and integration finished, pass
  else
    printf "\b== PASSED ==\n" | tee -a "${results_file}" | color $green
  fi

  # Empty temporary file
  > ${tmp_file}

  # Printing new line
  printf "\n" | tee -a "${results_file}"
}

# Counting number of tests
countries=0
total_tests=0
failed=0

# Iterating through folders in ..sample_data/
for folder in ../sample_data/*; do
  if [ -d "${folder}" ]; then
    countries=$((countries+1))
    country=${folder##*/}
    printf " -------- Testing ${country}\n\n" | tee -a "${results_file}" | color $magenta

  # Iterating through maps (GeoJSON(s) or JSON(s)) in country's folder
    for map in ${folder}/*.*json; do
      map_file_name=${map##*/}
      printf "Trying ${map_file_name}\n" | tee -a "${results_file}" | color $blue

      # Iterating through visual variable files (CSV(s)) in country's folder
      for csv in ${folder}/*.csv; do
        total_tests=$((total_tests+1))
        csv_name=${csv##*/}
        printf " - with ${csv_name}\n\n" | tee -a "${results_file}" | color $cyan

        # Running cartogram, and timing it in seconds
        run_map

      done
    done
    printf " -------- Finished testing ${country}.\n\n\n\n\n\n" | tee -a "${results_file}" | color $magenta
  fi
done

# Summary report
printf "===== Finished testing all countries. =====\n\n" | tee -a "${results_file}" | color $magenta

duration="$(($SECONDS / 60))m $(($SECONDS % 60))s"
printf "Finished ${total_tests} tests on ${countries} countries in ${duration}.\n"

failed_per=$(( 100 * $failed / $total_tests))
printf "Passed [$((total_tests-failed))/${total_tests}] | $((100-failed_per))%% \n" | tee -a "${results_file}" | color $green
printf "Failed [${failed}/${total_tests}] | ${failed_per}%% \n"  | tee -a "${results_file}" | color $red

# Checking if any tests failed
if [[ "${failed}" -gt 0 ]]; then

  # Showing failed tests and removing temporary file
  printf "\nFailed tests:\n"  | tee -a "${results_file}" | color $red
  cat failed_tmp.txt | tee -a "${results_file}"
  rm failed_tmp.txt
  printf "\n" | tee -a "${results_file}"
fi

# Removing temporary files
rm ${tmp_file}

# Prompting for file deletion
printf "\bClear ALL *.eps and *.geojson files in current directory? [y/N]: " | color $yellow
read to_clear
if [[ "$to_clear" == "y" ]]; then
  rm *.eps; rm *.geojson
  printf "All *.eps and *.geojson files deleted.\n" | color $red
else
  printf "Files not cleared.\n" | color $green
fi
