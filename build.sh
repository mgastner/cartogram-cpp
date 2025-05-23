# Find out if $BUILD_DIR is set, if not then assume it's in a local environment
if [[ -n ${BUILD_DIR} ]]; then
    # Docker is already a virtual environment so no need to create one
    # Build the project
    cmake -B "$BUILD_DIR"  -S . -DCMAKE_TOOLCHAIN_FILE=$BUILD_DIR/build/Release/generators/conan_toolchain.cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_EXPORT_COMPILE_COMMANDS=ON
    cmake --build "$BUILD_DIR" -j4
    cmake --install "$BUILD_DIR"
else
    # Create and activate virtual environment
    python3 -m venv .venv
    source .venv/bin/activate

    # Build the project
    sudo .venv/bin/cmake -B "build" -S . -DCMAKE_TOOLCHAIN_FILE=build/build/Release/generators/conan_toolchain.cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_EXPORT_COMPILE_COMMANDS=ON
    sudo .venv/bin/cmake --build "build" -j4
    sudo .venv/bin/cmake --install "build"
fi