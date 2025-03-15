# Install cartogram-cpp program
cmake -B build-docker
make -j1 -C build-docker
make install -C build-docker
echo "---PROGRAM BUILD DONE---"