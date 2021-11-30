#!/usr/bin/bash

WORK_DIR=$PWD

apt-get update && apt-get install -y --no-install-recommends \
        build-essential \
        libboost-all-dev \
        cmake \
        git

git clone https://github.com/andrewguy/fast-mRMR.git
cd fast-mRMR/cpu/src && make
cp fast-mrmr /usr/local/bin/fast-mrmr
cd ${WORK_DIR}/fast-mRMR/utils/data-reader && make
cp mrmr-reader /usr/local/bin/mrmr-reader
cd ${WORK_DIR}
rm -r ./fast-mRMR

git clone --recursive https://github.com/Jokeren/gBolt.git
cd gBolt && mkdir build && cd build && cmake .. && make
cp gbolt /usr/local/bin/gbolt
cd ${WORK_DIR}
rm -r ./gBolt