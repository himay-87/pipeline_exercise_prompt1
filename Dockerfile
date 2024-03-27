# Use a base image with necessary tools
FROM ubuntu:latest

# Install required packages
RUN apt-get update && apt-get install -y \
    bwa \
    samtools \
    picard-tools \
    freebayes \
    bcftools \
    perl \
    curl \
    unzip \
    default-jre \
    python3 \
    python3-pip \
    && rm -rf /var/lib/apt/lists/*

# Install VEP
RUN curl -LO https://github.com/Ensembl/ensembl-vep/archive/refs/tags/release/105.2.zip \
    && unzip 105.2.zip \
    && mv ensembl-vep-release-105.2 vep \
    && rm 105.2.zip \
    && cd vep \
    && perl INSTALL.pl --AUTO ap --PLUGINS Downstream,LoF,ExAC,SpliceRegion

# Set up environment variables for VEP
ENV PATH="/vep:${PATH}"
ENV PERL5LIB="/vep:${PERL5LIB}"

# Install FastQC and MultiQC
RUN apt-get update && apt-get install -y fastqc \
    && pip3 install multiqc

# Copy Picard Tools jar file
COPY picard.jar /picard.jar

# Set the working directory
WORKDIR /

# Copy the variant calling pipeline script and input files
COPY variant_calling_pipeline.py /variant_calling_pipeline.py
COPY reference.fasta /reference.fasta
COPY reads1.fastq /reads1.fastq
COPY reads2.fastq /reads2.fastq

# Documenting the Dockerfile

# This Dockerfile builds an image for a variant calling pipeline. It installs necessary tools and dependencies, 
# copies the pipeline script and input files, and sets up the environment for execution.

# - `FROM ubuntu:latest`: Uses the latest Ubuntu base image.
# - Installs required packages using `apt-get`.
# - Downloads and installs VEP (Variant Effect Predictor) for variant annotation.
# - Sets up environment variables for VEP.
# - Installs FastQC and MultiQC for quality control and analysis.
# - Copies the Picard Tools JAR file required for marking duplicates.
# - Sets the working directory to root.
# - Copies the variant calling pipeline script (`variant_calling_pipeline.py`) and input files (`reference.fasta`, `reads1.fastq`, `reads2.fastq`).
# - The final command (`CMD`) runs the variant calling pipeline script using Python 3.

# To build the Docker image, use:
#   docker build -t variant_calling_image .

# To run the Docker container using the built image, use:
#   docker run -it --rm variant_calling_image

# After running the container, the variant calling pipeline will execute automatically.
# CMD command to run the variant calling pipeline script
CMD ["python3", "variant_calling_pipeline.py"]
