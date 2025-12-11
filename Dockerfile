# Stage 1: Build dotplotrs from source (not available on conda)
FROM rust:slim AS dotplotrs-builder

RUN apt-get update && apt-get install -y \
    git \
    pkg-config \
    libfreetype6-dev \
    libfontconfig1-dev \
    && rm -rf /var/lib/apt/lists/*

RUN git clone https://github.com/Genoscope-LBGB/dotplotrs /tmp/dotplotrs \
    && cd /tmp/dotplotrs \
    && cargo build --release

# Stage 2: Final image with conda/mamba
FROM mambaorg/micromamba:1.5-jammy

USER root

# Copy dotplotrs binary from builder stage
COPY --from=dotplotrs-builder /tmp/dotplotrs/target/release/dotplotrs /usr/local/bin/

# Install runtime dependencies for dotplotrs
RUN apt-get update && apt-get install -y --no-install-recommends \
    libfreetype6 \
    libfontconfig1 \
    && rm -rf /var/lib/apt/lists/*

# Switch back to micromamba user
USER $MAMBA_USER

# Install conda packages
RUN micromamba install -y -n base -c conda-forge -c bioconda \
    python=3.11 \
    snakemake \
    biopython \
    xopen \
    taxonkit \
    mash \
    minimap2 \
    ncbi-datasets-cli \
    busco \
    && micromamba clean --all --yes

# Install Python packages not available on conda
ARG MAMBA_DOCKERFILE_ACTIVATE=1
RUN pip install --no-cache-dir \
    find_reference_genomes \
    snakemake-executor-plugin-slurm

# Install HoBRAC
USER root
WORKDIR /app
COPY . /app
RUN chown -R $MAMBA_USER:$MAMBA_USER /app

USER $MAMBA_USER
RUN pip install --no-cache-dir .

# Set environment variables
ENV TAXONKIT_DB=/taxonkit
