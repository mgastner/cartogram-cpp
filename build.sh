# Set the build directory from $BUILD_DIR if set, otherwise default to "build"
BUILD_DIR=${BUILD_DIR:-build}

cmake -B "$BUILD_DIR"
sudo make install -j4 -C "$BUILD_DIR"