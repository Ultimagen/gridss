ARG UBUNTU_VERSION=22.04
FROM ubuntu:$UBUNTU_VERSION AS gridss_base_closest_mirror
# Configure apt to use AWS EC2 mirrors for faster downloads on AWS infrastructure
RUN sed -i 's|http://archive.ubuntu.com/ubuntu/|http://us-east-1.ec2.archive.ubuntu.com/ubuntu/|g' /etc/apt/sources.list && \
    sed -i 's|http://security.ubuntu.com/ubuntu/|http://us-east-1.ec2.archive.ubuntu.com/ubuntu/|g' /etc/apt/sources.list

# Set up a C build environment for gridsstools, samtools, and R packages
FROM gridss_base_closest_mirror AS gridss_c_build_environment
RUN apt-get update && DEBIAN_FRONTEND=noninteractive apt-get install --no-install-recommends -y ca-certificates
RUN apt-get update && DEBIAN_FRONTEND=noninteractive apt-get install --no-install-recommends -y \
	libssl-dev \
	libcurl4-openssl-dev \
	libxml2-dev \
	zlib1g-dev \
	libbz2-dev \
	liblzma-dev \
	libdeflate-dev \
	build-essential \
	autotools-dev \
	autoconf \
	autogen \
	automake \
	libtool \
	make \
	wget \
	libomp-dev \
    python3-pip \
	&& rm -rf /var/lib/apt/lists/* \
	&& pip3 install pysam==0.22.0 biopython==1.83 pyfaidx==0.8.1.1 joblib==1.4.0 parasail==1.3.4 \
	&& ln -s /usr/bin/python3 /usr/bin/python

# compile gridsstools
FROM gridss_c_build_environment AS gridss_builder_c
RUN mkdir /opt/gridss/
COPY src/main/c /opt/gridss/src/main/c
COPY src/test/resources /opt/gridss/src/test/resources
RUN cd /opt/gridss/src/main/c/gridsstools/htslib && \
	autoreconf -i && ./configure && make -j 8 && \
	cd .. && \
	autoreconf -i && ./configure && make -j 8 && \
	cp gridsstools /opt/gridss/

# compile GRIDSS Java code
FROM maven:3.8.4-jdk-11 AS gridss_builder_java
RUN mkdir /opt/gridss/
WORKDIR /opt/gridss/
# Download maven dependencies first so docker can cache them
COPY pom.xml /opt/gridss/
COPY repo /opt/gridss/repo
# run all stages so all dependencies are cached
RUN mvn -Dmaven.artifact.threads=8 verify && rm -rf target
# Build GRIDSS jar
ARG GRIDSS_VERSION
COPY src /opt/gridss/src
RUN mvn -T 1C -Drevision=${GRIDSS_VERSION} package -Dmaven.test.skip=true && \
	cp target/gridss-${GRIDSS_VERSION}-gridss-jar-with-dependencies.jar /opt/gridss/

FROM gridss_c_build_environment AS gridss
# Setup CRAN ubuntu package repository
# apt-get clean not required for ubuntu images
# Pre-configure timezone to prevent interactive prompts
ENV TZ=UTC
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone
RUN apt-get update && DEBIAN_FRONTEND=noninteractive apt-get install --no-install-recommends -y \
		apt-transport-https \
		software-properties-common \
		dirmngr \
		gnupg && \
	apt-key adv \
	--keyserver hkp://keyserver.ubuntu.com:80 \
	--recv-keys 0xE298A3A825C0D65DFD57CBB651716619E084DAB9 && \
	add-apt-repository 'deb https://cloud.r-project.org/bin/linux/ubuntu jammy-cran40/' && \
	apt-get update && apt-get install --no-install-recommends -y \
		apt-utils \
		gawk \
		openjdk-17-jre-headless \
		bwa \
		hmmer \
		bedtools \
		bcftools \
		r-base \
		time \
		perl-modules \
		libtext-soundex-perl \
		python3-h5py \
		rsync \
		curl \
		libcairo2-dev \
		libgit2-dev \
		default-libmysqlclient-dev \
		libpq-dev \
		libsasl2-dev \
		libsqlite3-dev \
		libssh2-1-dev \
		libxtst6 \
		libharfbuzz-dev \
		libfribidi-dev \
		libfreetype6-dev \
		libpng-dev \
		libtiff5-dev \
		libjpeg-dev \
		libwebp-dev \
		libuv1-dev \
		unixodbc-dev \
        libncurses5-dev \
        libncursesw5-dev \
		gfortran \
		cmake \
		liblapack-dev \
		libblas-dev \
	&& rm -rf /var/lib/apt/lists/*

# samtools needs to be installed from source since the OS package verion is too old
RUN mkdir /opt/samtools && \
	cd /opt/samtools && \
	wget https://github.com/samtools/samtools/releases/download/1.14/samtools-1.14.tar.bz2 && \
	tar -jxf samtools-1.14.tar.bz2 && \
	cd samtools-1.14 && \
	autoheader && \
	autoconf -Wno-syntax && \
	./configure && \
	make install && \
	cd ~ && \
	rm -rf /opt/samtools
### Repeat Masker and dependencies
RUN mkdir /opt/trf && \
	cd /opt/trf && \
	wget https://github.com/Benson-Genomics-Lab/TRF/releases/download/v4.09.1/trf409.linux64 && \
	chmod +x trf*.linux64 && \
	ln -s trf*.linux64 trf
# Turns out we need makeblastdb as well as rmblastn (https://github.com/PapenfussLab/gridss/issues/535)
RUN mkdir /opt/rmblast
COPY rmblast-2.14.1+-x64-linux.tar.gz /opt/rmblast
RUN cd /opt/rmblast && \
	tar --no-anchored --strip-components 2 -xvzf rmblast-2.14.1+-x64-linux.tar.gz rmblastn makeblastdb && \
	rm rmblast-2.14.1+-x64-linux.tar.gz
COPY RepeatMasker-4.1.2-p1.tar.gz /opt/
RUN cd /opt/ && \
	tar zxf RepeatMasker-*.tar.gz && \
	rm RepeatMasker-*.tar.gz
# Install GATK
# GATK 4.7.0.0 rebuilt by Ultima with netty-bom 4.1.139.Final (CVE-2026-75595, no log4j 1.2.17), see BIOIN-3161.
# The spark jar is removed: only non-Spark tools (IntervalListToBed, PrintReads) are used.
RUN export GATK_VERSION="4.7.0.0-ug.netty1" && \
    export GATK_SHA256="044b98ad949d5a93eaf7666b5c2b5145741f68d7fd6a4fe75da109b83be02f78" && \
    mkdir -p /opt/gatk && \
    cd /opt/gatk && \
    wget https://github.com/Ultimagen/gatk/releases/download/${GATK_VERSION}/gatk-${GATK_VERSION}.zip && \
    echo "${GATK_SHA256}  gatk-${GATK_VERSION}.zip" | sha256sum -c - && \
    unzip gatk-${GATK_VERSION}.zip && \
    rm gatk-${GATK_VERSION}.zip && \
    rm gatk-${GATK_VERSION}/gatk-package-${GATK_VERSION}-spark.jar && \
    ln -s /opt/gatk/gatk-${GATK_VERSION}/gatk /usr/local/bin/gatk
### Kraken2 and dependencies
# dustmasker from e-direct: (or is this in ncbi-blast as well?)
RUN mkdir /opt/blast && \
	cd /opt/blast && \
	wget https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.11.0/ncbi-blast-2.11.0+-x64-linux.tar.gz && \
	tar zxf ncbi-blast-*.tar.gz && \
	mv ncbi-blast-*/bin/* . && \
	rm -r ncbi-blast-*
ENV KRAKEN_VERSION=2.1.2
RUN mkdir /opt/kraken2 && \
	cd /opt/kraken2 && \
	wget https://github.com/DerrickWood/kraken2/archive/refs/tags/v$KRAKEN_VERSION.tar.gz && \
	tar zxf v*.tar.gz && \
	cd kraken2* && \
	./install_kraken2.sh /opt/kraken2 && \
	cd .. && \
	rm -r kraken2-$KRAKEN_VERSION v*.tar.gz
RUN sh -c "$(wget -q ftp://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/install-edirect.sh -O -)" && \
	mv $HOME/edirect /opt/edirect
ENV PATH="/opt/gridss/:/opt/RepeatMasker:/opt/rmblast/:/opt/trf:/opt/kraken2:/opt/blast:/opt/edirect:/opt/gatk:$PATH"
# configure repeatmasker
RUN cd /opt/RepeatMasker && \
	perl configure \
		-default_search_engine rmblast \
		-rmblast_dir /opt/rmblast \
		-trf_prgm /opt/trf/trf \
		-hmmer_dir /usr/local/bin
# R packages used by GRIDSS - R package need the C toolchain installed
ENV R_INSTALL_STAGED=false
RUN Rscript -e 'options(Ncpus=8L, repos="https://cloud.r-project.org/", warn=2); \
	packages <- c("tidyverse", "assertthat", "testthat", "stringdist", "stringr", "argparser", "BiocManager", "pbapply"); \
	install.packages(packages); \
	failed <- packages[!packages %in% installed.packages()[,"Package"]]; \
	if(length(failed) > 0) { \
		cat("ERROR: Failed to install packages:", paste(failed, collapse=", "), "\n"); \
		quit(status=1); \
	}'
RUN Rscript -e 'options(Ncpus=8L, repos="https://cloud.r-project.org/", warn=2); \
	library(BiocManager); \
	packages <- c("StructuralVariantAnnotation", "VariantAnnotation", "rtracklayer", "BSgenome", "Rsamtools"); \
	BiocManager::install(ask=FALSE, pkgs=packages); \
	failed <- packages[!packages %in% installed.packages()[,"Package"]]; \
	if(length(failed) > 0) { \
		cat("ERROR: Failed to install Bioconductor packages:", paste(failed, collapse=", "), "\n"); \
		quit(status=1); \
	}'
# Install GRIDSS
ARG GRIDSS_VERSION
ENV GRIDSS_VERSION=${GRIDSS_VERSION}
ENV GRIDSS_JAR=/opt/gridss/gridss-${GRIDSS_VERSION}-gridss-jar-with-dependencies.jar
LABEL software="GRIDSS"
LABEL software.version="$GRIDSS_VERSION"
LABEL about.summary="Genomic Rearrangement IDentification Software Suite"
LABEL about.home="https://github.com/PapenfussLab/gridss"
LABEL about.tags="Genomics"
RUN mkdir /opt/gridss/ /data
COPY --from=gridss_builder_c /opt/gridss/gridsstools /opt/gridss/
COPY --from=gridss_builder_java /opt/gridss/gridss-${GRIDSS_VERSION}-gridss-jar-with-dependencies.jar /opt/gridss/
RUN pip3 install tqdm
# linux-libc-dev (kernel headers) only carries unfixed kernel CVEs and is needed only for compiling.
# Everything is compiled by now, so remove it together with the packages depending on it (libc6-dev,
# build-essential, *-dev). No autoremove: it would drop runtime libraries the R packages link against.
RUN apt-get purge -y linux-libc-dev && \
	rm -rf /var/lib/apt/lists/*
COPY scripts/gridss \
	scripts/gridss_annotate_vcf_kraken2 \
	scripts/gridss_annotate_vcf_repeatmasker \
	scripts/gridss_extract_overlapping_fragments \
	scripts/gridss_somatic_filter \
	scripts/virusbreakend \
	scripts/virusbreakend-build \
	scripts/gridss.config.R \
	scripts/libgridss.R \
	scripts/link_breakpoints.R \
	scripts/choose_best_haplotype_realignment.py \
	scripts/align_long_homopolymers.py \
	scripts/rematch_reads_to_haplotypes.py \
	scripts/realigned_regions.py \
	scripts/convert_vcf_format.R \
	/opt/gridss/
RUN chmod +x /opt/gridss/* && \
	chmod -x /opt/gridss/*.R
WORKDIR /data/

# Include licenses
COPY ["LICENSE", "."]
COPY ["third_party_licenses/gatk LICENSE.txt", "."]


# Copy build artifact locally
FROM scratch AS gridss_export_build_artefacts
ARG GRIDSS_VERSION
COPY --from=gridss /opt/gridss/* /

