FROM python:3.12-slim-bookworm

# Install dependencies
RUN apt-get update
RUN apt install -y g++-11 build-essential cmake libboost-all-dev nlohmann-json3-dev libomp-dev libfftw3-dev libcairo2-dev libmpfr-dev libgmp-dev libboost-dev

# Set environment variable for build shell script to indicate inside Docker environment
ENV IN_DOCKER_CONTAINER=true

# Create and set working directory to "cartogram"
WORKDIR /cartogram