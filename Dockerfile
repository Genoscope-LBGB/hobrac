FROM python:3.11-slim

# Set non-interactive frontend
ENV DEBIAN_FRONTEND=noninteractive

# 1. Install system dependencies
# build-essential: for compiling tools
# curl/wget: for downloading binaries
# zlib1g-dev: often needed by bioinformatics tools
# git: for pip install git+... if needed
# pkg-config, libfreetype6-dev, libfontconfig1-dev: potential deps for plotting libs in Rust
RUN apt-get update && apt-get install -y \
    wget \
    curl \
    build-essential \
    zlib1g-dev \
    git \
    pkg-config \
    libfreetype6-dev \
    libfontconfig1-dev \
    unzip \
    && rm -rf /var/lib/apt/lists/*

# 2. Install Rust (for dotplotrs)
RUN curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh -s -- -y
ENV PATH="/root/.cargo/bin:${PATH}"

# 3. Install dotplotrs (clone and build from source)
RUN git clone https://github.com/Genoscope-LBGB/dotplotrs /tmp/dotplotrs \
    && cd /tmp/dotplotrs \
    && cargo build --release \
    && mv /tmp/dotplotrs/target/release/dotplotrs /usr/local/bin/ \
    && rm -rf /tmp/dotplotrs

# 4. Install Bioinformatics Tools (Binary Downloads)

# Taxonkit
ENV TAXONKIT_VERSION=v0.16.0
RUN wget https://github.com/shenwei356/taxonkit/releases/download/${TAXONKIT_VERSION}/taxonkit_linux_amd64.tar.gz \
    && tar -zxvf taxonkit_linux_amd64.tar.gz \
    && mv taxonkit /usr/local/bin/ \
    && rm taxonkit_linux_amd64.tar.gz

# Mash
ENV MASH_VERSION=v2.3
RUN wget https://github.com/marbl/Mash/releases/download/${MASH_VERSION}/mash-Linux64-${MASH_VERSION}.tar \
    && tar -xf mash-Linux64-${MASH_VERSION}.tar \
    && mv mash-Linux64-${MASH_VERSION}/mash /usr/local/bin/ \
    && rm -rf mash-Linux64-${MASH_VERSION}*

# Minimap2
ENV MINIMAP2_VERSION=2.30
RUN wget https://github.com/lh3/minimap2/releases/download/v${MINIMAP2_VERSION}/minimap2-${MINIMAP2_VERSION}_x64-linux.tar.bz2 \
    && tar -jxvf minimap2-${MINIMAP2_VERSION}_x64-linux.tar.bz2 \
    && mv minimap2-${MINIMAP2_VERSION}_x64-linux/minimap2 /usr/local/bin/ \
    && rm -rf minimap2-${MINIMAP2_VERSION}_x64-linux*

# NCBI Datasets & find_reference_genomes
# Installing find_reference_genomes via pip later, but datasets cli needs binary
RUN curl -o datasets 'https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/v2/linux-amd64/datasets' \
    && chmod +x datasets \
    && mv datasets /usr/local/bin/

# 5. Install Python Packages
RUN pip install --no-cache-dir \
    snakemake \
    biopython \
    xopen \
    find_reference_genomes

# 6. Install HoBRAC
# Copy the current directory contents into the container at /app
WORKDIR /app
COPY . /app
RUN pip install .

# 7. Set Environment Variables
ENV TAXONKIT_DB=/taxonkit

