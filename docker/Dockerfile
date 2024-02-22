FROM ubuntu:22.04

# Install dependencies
RUN apt-get update && apt-get install -y \
	python3.10 python3-pip \
    build-essential \
    manpages-dev \
    software-properties-common \
    libcgal-dev \
    libomp-dev \
    libcairo2-dev \
    nlohmann-json3-dev \
    fftw3-dev \
    && rm -rf /var/lib/apt/lists/*
	
# Update alternatives to use python3 as the default python version
RUN update-alternatives --install /usr/bin/python python /usr/bin/python3.10 1 \
    && update-alternatives --install /usr/bin/pip pip /usr/bin/pip3 1

# Copy the compiled binary
ARG CARTOGRAM_BIN
COPY ${CARTOGRAM_BIN} /opt/cartogram

CMD ["tail", "-f", "/dev/null"]