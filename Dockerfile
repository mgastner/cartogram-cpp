FROM python:3.12-slim-bookworm

# Set environment variable for build shell script to indicate inside Docker environment
ENV BUILD_DIR=build-docker

# Create and set working directory to "cartogram"
WORKDIR /cartogram

# Copy over all project files
COPY . .

# Remove Windows-style line endings in shell scripts
RUN sed -i 's/\r$//' build.sh
RUN sed -i 's/\r$//' tests/stress_test.sh

# Install Clang
RUN apt update
RUN apt install -y build-essential clang

# Install dependencies
RUN pip install --upgrade pip wheel conan==2.16.1 cmake==3.30.0

# Setup Conan
RUN conan remote update conancenter --url=https://center2.conan.io
RUN conan profile detect

# Install dependencies via Conan
RUN conan install . --output-folder build-docker --build=missing -s build_type=Release -s compiler.cppstd=20

# Build the project
RUN cmake -B build-docker  -S . -DCMAKE_TOOLCHAIN_FILE=build-docker/build/Release/generators/conan_toolchain.cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_EXPORT_COMPILE_COMMANDS=ON
RUN cmake --build build-docker -j4
RUN cmake --install build-docker

# Change working directory to output folder
WORKDIR /cartogram/output