# 使用Ubuntu 20.04基础镜像
FROM ubuntu:20.04

# 设置非交互式环境
ENV DEBIAN_FRONTEND=noninteractive

# 安装系统工具（包括 wget 和 CA 证书）
RUN apt-get update && apt-get install -y \
    wget \
    unzip \
    default-jre \
    bzip2 \
    ca-certificates \
    && update-ca-certificates \
    && rm -rf /var/lib/apt/lists/*

# 安装 Miniconda
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda.sh && \
    bash ~/miniconda.sh -b -p /opt/conda && \
    rm ~/miniconda.sh

# 将 Conda 添加到 PATH
ENV PATH="/opt/conda/bin:$PATH"

# 配置 Conda 国内镜像源（此时 conda 已可用）
RUN conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/main && \
    conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/r && \
    conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/msys2 && \
    conda config --set show_channel_urls yes

# 安装其他依赖（通过 apt 或 conda）
RUN apt-get update && apt-get install -y \
    perl \
    python3 \
    python3-pip \
    build-essential \
    zlib1g-dev \
    libncurses5-dev \
    libbz2-dev \
    liblzma-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    && rm -rf /var/lib/apt/lists/*

# 使用 mamba 安装生信工具（指定 Python 3.11 和工具版本）
RUN conda install -y -n base python=3.11 && \
    conda install -y -c conda-forge -c bioconda \
    snakemake \
    sra-tools=3.0.10 \
    fastqc=0.12.1 \
    trimmomatic=0.39 \
    bwa=0.7.17 \
    samtools=1.16.1 \
    bcftools=1.16.0 && \
    conda clean -a -y

# 设置Trimmomatic环境变量（conda安装的路径）
ENV TRIMMOMATIC_JAR=/opt/conda/share/trimmomatic/trimmomatic.jar

# 创建Trimmomatic快捷方式
RUN echo '#!/bin/bash\njava -jar $TRIMMOMATIC_JAR "$@"' > /usr/local/bin/trimmomatic && \
    chmod +x /usr/local/bin/trimmomatic

# 确保系统时间正确（解决SSL证书验证问题）
RUN apt-get update && apt-get install -y tzdata && \
    ln -fs /usr/share/zoneinfo/UTC /etc/localtime && \
    dpkg-reconfigure --frontend noninteractive tzdata && \
    rm -rf /var/lib/apt/lists/*
