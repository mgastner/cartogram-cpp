cmake .
make
cd bin
./cartogram -g ../sample_data/belgium/belgium.geojson ../sample_data/belgium/belgium_population2019.csv -e
cd ..