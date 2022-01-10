cmake .
make
./bin/cartogram -g ./sample_data/belgium/belgium.geojson ./sample_data/belgium/belgium_population2019.csv -e
./bin/cartogram -g ./sample_data/vietnam/vietnam.geojson ./sample_data/vietnam/vietnam_population2019.csv -e