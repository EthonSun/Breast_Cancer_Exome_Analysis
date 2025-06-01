

# 乳腺癌外显子测序分析流程 (Docker & Snakemake)

这是一个利用 **Docker** 和 **Snakemake** 实现的乳腺癌外显子测序数据分析完整流程，旨在提供一个可重复、自动化且易于部署的生物信息学分析环境。

---

## 1. 准备工作目录结构

```bash
# 创建主项目目录
mkdir -p ~/course/2503/project
bio=${HOME}/course/2503/project

# 创建多层目录结构
mkdir -p ${bio}/{data,metadata,logs,results,scripts}
mkdir -p ${bio}/data/{sra,ref,trimmed_fastq,untrimmed_fastq}
mkdir -p ${bio}/results/{fastqc_untrimmed,fastqc_trimmed,sam,bam,bcf,vcf}
```

------

## 2. 构建Docker镜像

首先，创建用于存放 Dockerfile 的目录，并进入该目录：

```bash
# 创建Dockerfile目录
mkdir ~/linux_docker
cd ~/linux_docker
```

### 创建Dockerfile文件：

将以下内容保存为 `Dockerfile`：

```dockerfile
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
```

### 构建镜像：

在 `~/linux_docker` 目录下执行以下命令构建 Docker 镜像：

```bash
cd /path/to/your/project/docker 
docker build -t linux_ngs_pipeline:latest .

# 或者直接在Docker Hub中拉取
docker pull ethonsun/breast_cancer_pipeline:latest
```

------

## 3. 下载参考基因组

在宿主机上下载参考基因组并创建索引：

```bash
# 在宿主机上下载参考基因组
cd ${bio}/data/ref

# 下载hg38参考基因组
wget [https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz](https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz)
gunzip hg38.fa.gz
mv hg38.fa Homo_sapiens_assembly38.fasta

# 创建索引（在Docker容器内完成）
docker run -it --rm -v ${bio}:/data linux_ngs_pipeline:latest \
    bash -c "cd /data/data/ref && \
    bwa index Homo_sapiens_assembly38.fasta && \
    samtools faidx Homo_sapiens_assembly38.fasta"
```

------

## 4. 下载SRA数据

创建下载脚本 `${bio}/scripts/download_sra.sh`，并复制以下内容：

```bash
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
```

运行下载：

```bash
docker run -it --rm -v ${bio}:/data linux_ngs_pipeline:latest \
    bash /data/scripts/download_sra.sh
```

------

## 5. 创建分析流程脚本 (Snakemake)

创建 `${bio}/scripts/Snakefile` 文件，并复制以下内容：

```python
import os
import glob

configfile: "config.yaml"

# 获取样本列表（使用os.path兼容路径）
SAMPLES = [os.path.basename(f).replace('.sra','') 
          for f in glob.glob(os.path.join(config['sra_dir'], "*", "*.sra"))]

rule all:
    input:
        expand(os.path.join(config['result_dir'], "vcf/{sample}_final.vcf"), sample=SAMPLES)

# 1. SRA转fastq
rule sra_to_fastq:
    input:
        sra = os.path.join(config['sra_dir'], "{sample}", "{sample}.sra")
    output:
        r1 = os.path.join(config['untrimmed_dir'], "{sample}_1.fastq.gz"),
        r2 = os.path.join(config['untrimmed_dir'], "{sample}_2.fastq.gz")
    log:
        log = os.path.join("logs/sra_to_fastq", "{sample}.log")
    params:
        outdir = config['untrimmed_dir']
    shell:
        """
        mkdir -p {params.outdir}
        fasterq-dump {input.sra} -O {params.outdir} &> {log.log}
        gzip {params.outdir}/{wildcards.sample}_1.fastq
        gzip {params.outdir}/{wildcards.sample}_2.fastq
        """

# 2. FastQC原始数据
rule fastqc_raw:
    input:
        r1 = rules.sra_to_fastq.output.r1,
        r2 = rules.sra_to_fastq.output.r2
    output:
        html1 = os.path.join(config['fastqc_raw_dir'], "{sample}_1_fastqc.html"),
        html2 = os.path.join(config['fastqc_raw_dir'], "{sample}_2_fastqc.html"),
        zip1 = os.path.join(config['fastqc_raw_dir'], "{sample}_1_fastqc.zip"),
        zip2 = os.path.join(config['fastqc_raw_dir'], "{sample}_2_fastqc.zip")
    log:
        log = os.path.join("logs/fastqc_raw", "{sample}.log")
    params:
        outdir = config['fastqc_raw_dir']
    shell:
        """
        mkdir -p {params.outdir} logs/fastqc_raw
        fastqc {input.r1} {input.r2} -o {params.outdir} &> {log.log}
        """

# 3. Trimming
rule trim:
    input:
        r1 = rules.fastqc_raw.input.r1,
        r2 = rules.fastqc_raw.input.r2
    output:
        r1_paired = os.path.join(config['trimmed_dir'], "{sample}_1_paired.fastq.gz"),
        r1_unpaired = os.path.join(config['trimmed_dir'], "{sample}_1_unpaired.fastq.gz"),
        r2_paired = os.path.join(config['trimmed_dir'], "{sample}_2_paired.fastq.gz"),
        r2_unpaired = os.path.join(config['trimmed_dir'], "{sample}_2_unpaired.fastq.gz")
    log:
        log = os.path.join("logs/trim", "{sample}.log")
    params:
        outdir = config['trimmed_dir']
    threads: 2
    shell:
        """
        mkdir -p {params.outdir}
        trimmomatic PE -threads {threads} {input.r1} {input.r2} \
        {output.r1_paired} {output.r1_unpaired} \
        {output.r2_paired} {output.r2_unpaired} \
        SLIDINGWINDOW:4:20 MINLEN:20 &> {log.log}
        """

# 4. BWA比对
rule bwa_mem:
    input:
        r1 = rules.trim.output.r1_paired,
        r2 = rules.trim.output.r2_paired,
        ref = config['ref']
    output:
        sam = os.path.join(config['result_dir'], "sam/{sample}.sam")
    log:
        main = os.path.join("logs", "{sample}.bwa.log")
    threads: 2
    shell:
        """
        bwa mem -t {threads} -R '@RG\\tID:{wildcards.sample}\\tSM:{wildcards.sample}\\tPL:Illumina' \
        {input.ref} {input.r1} {input.r2} > {output.sam} 2> {log.main}
        """

# 5. SAM转BAM
rule sam_to_bam:
    input:
        rules.bwa_mem.output.sam
    output:
        bam = os.path.join(config['result_dir'], "bam/{sample}.sorted.bam"),
        bai = os.path.join(config['result_dir'], "bam/{sample}.sorted.bam.bai")
    log:
        main = os.path.join("logs", "{sample}.samtools.log")
    threads: 2
    shell:
        """
        samtools view -S -b {input} 2>> {log.main} | \
        samtools sort -@ {threads} -o {output.bam} 2>> {log.main} && \
        samtools index {output.bam} 2>> {log.main}
        """

# 6. 变异检测
rule call_variants:
    input:
        bam = rules.sam_to_bam.output.bam,
        ref = config['ref']
    output:
        vcf = os.path.join(config['result_dir'], "vcf/{sample}_final.vcf")
    log:
        main = os.path.join("logs", "{sample}.bcftools.log")
    shell:
        """
        bcftools mpileup -O b -f {input.ref} {input.bam} 2>> {log.main} | \
        bcftools call --ploidy 1 -m -v - 2>> {log.main} | \
        vcfutils.pl varFilter > {output.vcf} 2>> {log.main}
        """
```

配套配置文件 `${bio}/scripts/config.yaml`：

```yaml
sra_dir: "/data/data/sra"
untrimmed_dir: "/data/data/untrimmed"
trimmed_dir: "/data/data/trimmed"
ref: "/data/data/ref/Homo_sapiens_assembly38.fasta"
result_dir: "/data/results"
fastqc_raw_dir: "/data/results/fastqc_untrimmed"
fastqc_trimmed_dir: "/data/results/fastqc_trimmed"
```

### 执行流程：

在 Docker 容器内执行 Snakemake 流程：

```bash
docker run -it --rm \
  -v ~/course/2503/project:/data \
  linux_ngs_pipeline:latest \
  /bin/bash
  
# 进入容器后执行
cd /data/scripts
snakemake -j 4 --use-conda # 生成最终结果.vcf文件
```

------

## 6. 结果目录结构

流程完成后，`project/` 目录将包含以下结构和文件：

```
project/
├── data/
│   ├── Homo_sapiens_assembly38.fasta
│   ├── Homo_sapiens_assembly38.fasta.amb
│   ├── Homo_sapiens_assembly38.fasta.ann
│   ├── Homo_sapiens_assembly38.fasta.bwt
│   ├── Homo_sapiens_assembly38.fasta.fai
│   ├── Homo_sapiens_assembly38.fasta.pac
│   ├── Homo_sapiens_assembly38.fasta.sa
│   ├── sra/
│   │   └── SRR21389034/
│   │       └── SRR21389034.sra
│   ├── trimmed/
│   │   ├── SRR21389034_1_paired.fastq.gz
│   │   ├── SRR21389034_1_unpaired.fastq.gz
│   │   ├── SRR21389034_2_paired.fastq.gz
│   │   └── SRR21389034_2_unpaired.fastq.gz
│   ├── trimmed_fastq/
│   ├── untrimmed/
│   │   ├── SRR21389034_1.fastq.gz
│   │   └── SRR21389034_2.fastq.gz
│   └── untrimmed_fastq/
├── download.log
├── logs/
│   ├── SRR21389034.bwamen.log
│   ├── SRR21389034.log
│   ├── pipeline_20250531_215903.log
│   └── pipeline_20250601_110058.log
├── metadata/
├── results/
│   ├── bam/
│   │   ├── SRR21389034.bam
│   │   ├── SRR21389034.sorted.bam
│   │   └── SRR21389034.sorted.bam.bai
│   ├── bcf/
│   ├── fastqc_trimmed/
│   │   ├── SRR21389034_1_paired_fastqc.html
│   │   ├── SRR21389034_1_paired_fastqc.zip
│   │   ├── SRR21389034_2_paired_fastqc.html
│   │   └── SRR21389034_2_paired_fastqc.zip
│   ├── fastqc_untrimmed/
│   │   ├── SRR21389034_1_fastqc.html
│   │   ├── SRR21389034_1_fastqc.zip
│   │   ├── SRR21389034_2_fastqc.html
│   │   └── SRR21389034_2_fastqc.zip
│   ├── sam/
│   │   └── SRR21389034.sam
│   └── vcf/
│       └── SRR21389034_final.vcf
└── scripts/
    ├── Snakefile
    ├── config.yaml
    ├── download_sra.sh
    └── logs/
        ├── SRR21389034.bcftools.log
        └── SRR21389034.samtools.log
```

------

## 关键要点说明

- **Docker环境：** 完美解决了环境依赖问题，确保所有工具正常运行。
- **参考基因组：** 使用 UCSC hg38 版本，包含所有必要索引文件。
- 流程优化：
  - 使用多线程加速处理（根据系统资源调整）。
  - 每个步骤均记录详细日志。
  - 自动完成从原始数据到变异检测的全流程。

------

## 常见问题处理

### 磁盘空间不足：

- 清理中间文件：

  Bash

  ```
  rm ${bio}/data/sra/SRR*/SRR*.sra
  rm ${bio}/results/sam/*.sam
  ```

### 流程中断后继续：

- **删除未完成的输出文件后重新提交。**

### 查看运行状态：

- `qstat -u $USER` (如果使用 PBS/Slurm 作业系统)
- `docker ps -a`

------

## 许可证 (可选)

本项目采用 [MIT License](https://opensource.org/licenses/MIT) 。

------

## 作者

EthonSun/孙义宸

zihdeng/邓梓衡
