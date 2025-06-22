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

if [[ -z ${BUILD_TYPE} ]]; then
    BUILD_TYPE=Release
fi

# Output BUILD_DIR and BUILD_TYPE
echo "Building in directory: $BUILD_DIR, with build type: $BUILD_TYPE"

# Warn user if BUILD_DIR already exists
if [[ -d $BUILD_DIR ]]; then
    echo "Warning: $BUILD_DIR already exists. Pre-existing cache may interfere with the build process and cause installation to fail. Consider deleting it before proceeding."
fi

conan install . --output-folder $BUILD_DIR --build=missing -s build_type=$BUILD_TYPE -s compiler.cppstd=20

cmake -B "$BUILD_DIR"  -S . -DCMAKE_TOOLCHAIN_FILE=$BUILD_DIR/build/$BUILD_TYPE/generators/conan_toolchain.cmake -DCMAKE_BUILD_TYPE=$BUILD_TYPE -DCMAKE_EXPORT_COMPILE_COMMANDS=ON
cmake --build "$BUILD_DIR" -j4
cmake --install "$BUILD_DIR"
