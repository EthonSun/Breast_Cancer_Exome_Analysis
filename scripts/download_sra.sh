#!/bin/bash
# download_sra.sh

WORKDIR=/data
SRA_DIR=${WORKDIR}/data/sra
mkdir -p ${SRA_DIR}

# 添加--progress和--verbose参数
prefetch SRR21389034 \
    --progress \
    --verbose \
    -O ${SRA_DIR} 2>&1 | tee ${WORKDIR}/download.log
