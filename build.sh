# Check if inside Docker or non-Docker environment
if [[ "$IN_DOCKER_CONTAINER" && "$IN_DOCKER_CONTAINER" == "true" ]]; then
    # Build and install program using build-docker folder if set inside Docker container
    cmake -B build-docker
    make -j1 -C build-docker
    make install -C build-docker
else
    # Build and install program using build folder if not set inside Docker container
    cmake -B build
    make -C build
    sudo make install -C build
fi

echo "---PROGRAM BUILD DONE---"