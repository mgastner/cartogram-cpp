# Set BUILD_DIR and BUILD_TYPE if empty
if [[ -z ${BUILD_DIR} ]]; then
    BUILD_DIR=build
fi

if [[ -z ${BUILD_TYPE} ]]; then
    BUILD_TYPE=Release
fi

# Output BUILD_DIR and BUILD_TYPE
echo "Building in directory: $BUILD_DIR, with build type: $BUILD_TYPE"

# Warn user if BUILD_DIR already exists
if [[ -d $BUILD_DIR ]]; then
    echo "Warning: $BUILD_DIR already exists. Pre-existing cache may interfere with the build process and cause installation to fail. Consider deleting it before proceeding."
fi

pipx run conan==2.16.1 install . --output-folder $BUILD_DIR --build=missing -s build_type=$BUILD_TYPE -s compiler.cppstd=20

pipx run cmake==3.30.0 -B "$BUILD_DIR"  -S . -DCMAKE_TOOLCHAIN_FILE=$BUILD_DIR/build/$BUILD_TYPE/generators/conan_toolchain.cmake -DCMAKE_BUILD_TYPE=$BUILD_TYPE -DCMAKE_EXPORT_COMPILE_COMMANDS=ON
pipx run cmake==3.30.0 --build "$BUILD_DIR" -j4
pipx run cmake==3.30.0 --install "$BUILD_DIR"