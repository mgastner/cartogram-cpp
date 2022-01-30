#!/usr/bin/env bash

# Turning on extended glob option
shopt -s extglob

# Start time, and other metadata
start_date=$(date '+%Y-%m-%d_%H-%M-%S')
results="results_${start_date}.txt"
SECONDS=0
echo -e "Writing to ${results}\n"
echo -e "Tested on ${start_date}\n" >> "${results}"

# Parsing command line arguments
# From https://unix.stackexchange.com/questions/462515/
# check-command-line-arguments
if [ $# -eq 0 ] || [ -z "$1" ]; then
  cli="-est"
else
  cli="$1"
fi
echo -e "Passing ${cli} to cartogram\n" | tee -a "${results}"

# Counting number of tests
countries=0
total_tests=0
failed=0

# Iterating through folders in ..sample_data/
for folder in ../sample_data/*; do
  if [ -d "${folder}" ]; then
    countries=$((countries+1))
    country=${folder##*/}
    echo -e "-------- Testing ${country}...\n" | tee -a "${results}"

  # Iterating through maps (GeoJSON(s) or JSON(s)) in country's folder
    for map in ${folder}/*.*json; do
      map_file_name=${map##*/}
      echo "Trying ${map_file_name} ..." | tee -a "${results}"

      # Iterating through visual variable files (CSV(s)) in country's folder
      for csv in ${folder}/*.csv; do
        total_tests=$((total_tests+1))
        csv_name=${csv##*/}
        echo -e "with ${csv_name} ...\n" | tee -a "${results}"

        # Running cartogram, and timing it in seconds
        start=$(date +%s)
        $(cartogram ${map} ${csv} ${cli} > tmp.txt 2>&1) >> tmp.txt
        end=$(date +%s)
        runtime=$((end-start))

        # Checking for errors
        grep -i -A 3 "error" tmp.txt | tee -a "${results}"
        grep -i "warning" tmp.txt | tee -a "${results}"
        grep -i -m 1 -A 1 "invalid" tmp.txt | tee -a "${results}"
        echo -e "\nRuntime: ${runtime}s" | tee -a "${results}"

        # Checking if integration finished
        if grep -Fxq "Progress: 1" tmp.txt; then
          echo "== Integration finished ==" | tee -a "${results}"
        fi

        # Checking for any errors, invalid geometry or unfinished integration
        if grep -qi "invalid" tmp.txt || grep -qi "error" tmp.txt || ! grep -Fxq "Progress: 1" tmp.txt ; then
          echo -e "== FAILED ==" | tee -a "${results}"

          # Printing country to failed_tmp.txt
          if [[ "${failed}" -eq 0 ]] || ! grep -qi "${country}" failed_tmp.txt; then
            echo -e "\n${country}" >> failed_tmp.txt
          fi

          # Increasing failed counter
          failed=$((failed+1))

          # Saving data to file if FAILED
          map_wo_ext=${map_file_name%.*json}
          csv_wo_ext=${csv_name%.csv}
          err_file="results_${start_date}-${map_wo_ext}-${csv_wo_ext}.txt"
          echo "- ${map_file_name} with ${csv_name}" >> failed_tmp.txt
          echo -e "Full output saved to ${err_file}" | tee -a "${results}"
          cp tmp.txt "${err_file}"

        # If no errors, and integration finished, pass
        else
          echo -e "== PASSED ==" | tee -a "${results}"
        fi

        # Printing new line
        echo "" | tee -a "${results}"
      done
    done
    echo -e "-------- Finished testing ${country}.\n\n\n\n\n" | tee -a "${results}"
  fi
done

# Summary report
echo -e "===== Finished testing all countries. =====\n"  | tee -a "${results}"
duration="$(($SECONDS / 60))m $(($SECONDS % 60))s"
echo "Finished ${total_tests} tests on ${countries} countries in ${duration}. Passed $((total_tests-failed)), failed ${failed}."  | tee -a "${results}"

# Checking if any tests failed
if [[ "${failed}" -gt 0 ]]; then

  # Showing failed tests and removing temporary file
  echo "Failed tests:"  | tee -a "${results}"
  cat failed_tmp.txt | tee -a "${results}"
  rm failed_tmp.txt
  echo "" | tee -a "${results}"
fi

# Removing temporary files
rm tmp.txt

# Prompting for file deletion
read -p "Clear ALL *.eps and *.geojson files in current directory? [y/N]: " to_clear
if [[ "$to_clear" == "y" ]]; then
  rm *.eps; rm *.geojson
  echo "All *.eps and *.geojson files deleted."
else
  echo "Files not cleared."
fi
