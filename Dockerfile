FROM python:3.12-slim-bookworm

# Install dependencies
RUN apt-get update
RUN apt install -y g++-11 build-essential cmake libboost-all-dev nlohmann-json3-dev libomp-dev libfftw3-dev libcairo2-dev libmpfr-dev libgmp-dev libboost-dev

# Create and set working directory to "cartogram"
WORKDIR /cartogram
# Copy over all project files
COPY . .

# Install cartogram-cpp program
RUN cmake -B build-docker
RUN make -C build-docker
RUN make install -C build-docker