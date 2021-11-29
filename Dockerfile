FROM ubuntu:20.04
LABEL maintainer="Andrew Guy"

ENV DEBIAN_FRONTEND noninteractive

WORKDIR /app

COPY . .

RUN apt-get update && apt-get install -y --no-install-recommends \
        python3-pip \
        build-essential \
        libboost-all-dev \
        cmake \
        git && \
    rm -rf /var/lib/apt/lists/* && \
    git clone https://github.com/andrewguy/fast-mRMR.git && \
    cd fast-mRMR/cpu/src && \
    make && \
    cp fast-mrmr /usr/local/bin/fast-mrmr && \
    cd ../../utils/data-reader && \
    make && \
    cp mrmr-reader /usr/local/bin/mrmr-reader && \
    cd ../../.. && \
    rm -r fast-mRMR && \
    git clone --recursive https://github.com/Jokeren/gBolt.git && \
    cd gBolt && \
    mkdir build && cd build && \
    cmake .. && \
    make && \
    cp gbolt /usr/local/bin/gbolt && \
    cd ../.. && rm -r gBolt && \
    pip3 install numpy \
        pandas \
        sklearn \
        matplotlib \
        jupyter \
        networkx \
        statsmodels \
        pyxdameraulevenshtein \
        click \
        pyparsing && \
    python3 /app/setup.py install && \
    apt-get purge -y --auto-remove build-essential cmake git && \
    rm -r /app

RUN useradd -l -m -s /bin/bash -N -u 1000 user && \
    mkdir /data && \
    chown 1000:100 /data

USER 1000

WORKDIR /data

ENTRYPOINT [ "ccarl" ]
CMD ["--help"]