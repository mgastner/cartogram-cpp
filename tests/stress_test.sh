#!/usr/bin/env bash

shopt -s extglob

for folder in ../sample_data/*; do
  if [ -d "${folder}" ]; then
    echo "Working on ${country}..."
    for map in ${folder}/*.*json; do
      echo ${map}
      for csv in ${folder}/*.csv; do
        echo ${csv}
        $(cartogram ${map} -V ${csv} -est)
      done
    done
  fi
  rm *.eps; rm *.geojson
  exit 1
done
