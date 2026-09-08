FROM r-base:4.4.2@sha256:fe9b29520eeb5292d814b0958783c0ddfcdab37402967a3e67307604354f98d7 AS r-packages

RUN apt-get update \
    && apt-get install --yes --no-install-recommends \
        build-essential \
        ca-certificates \
        libcurl4-openssl-dev \
        libssl-dev \
        libxml2-dev \
        zlib1g-dev \
    && rm -rf /var/lib/apt/lists/*

COPY install_r_dependencies.R /tmp/install_r_dependencies.R
RUN Rscript /tmp/install_r_dependencies.R \
    && rm -f /tmp/install_r_dependencies.R

FROM r-base:4.4.2@sha256:fe9b29520eeb5292d814b0958783c0ddfcdab37402967a3e67307604354f98d7 AS blast-tools

ARG DEBIAN_FRONTEND=noninteractive
ARG BLAST_VERSION=2.16.0
ARG BLAST_MD5=48f66c9e01ea5136e381b2bf6fc62036

RUN apt-get update \
    && apt-get install --yes --no-install-recommends ca-certificates curl \
    && rm -rf /var/lib/apt/lists/*

RUN mkdir -p /tmp/blast /opt/blast \
    && curl --fail --silent --show-error --location --retry 3 \
        "https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/${BLAST_VERSION}/ncbi-blast-${BLAST_VERSION}+-x64-linux.tar.gz" \
        --output /tmp/ncbi-blast.tar.gz \
    && printf '%s  %s\n' "${BLAST_MD5}" /tmp/ncbi-blast.tar.gz | md5sum --check - \
    && tar --extract --gzip --file /tmp/ncbi-blast.tar.gz --directory /tmp/blast --strip-components=1 \
    && install --mode=0755 /tmp/blast/bin/blastn /opt/blast/blastn \
    && install --mode=0755 /tmp/blast/bin/blastp /opt/blast/blastp \
    && install --mode=0755 /tmp/blast/bin/makeblastdb /opt/blast/makeblastdb

FROM r-base:4.4.2@sha256:fe9b29520eeb5292d814b0958783c0ddfcdab37402967a3e67307604354f98d7

ARG IMAGE_VERSION=0.2.0

LABEL org.opencontainers.image.title="Standalone pneumococcal beta-lactam MIC predictor" \
      org.opencontainers.image.version="${IMAGE_VERSION}" \
      org.opencontainers.image.source="https://github.com/pathogenwatch/spn-resistance-pbp" \
      org.opencontainers.image.description="Assembly FASTA PBP1A/PBP2B/PBP2X Random Forest predictor"

RUN apt-get update \
    && apt-get install --yes --no-install-recommends \
        clustalo \
        libjson-perl \
    && rm -rf /var/lib/apt/lists/*

COPY --from=r-packages /usr/local/lib/R/site-library/ /usr/local/lib/R/site-library/
COPY --from=blast-tools /opt/blast/ /opt/blast/

ENV PATH="/opt/blast:/predictor:${PATH}"

WORKDIR /predictor
COPY SPN_Reference_DB/ /predictor/SPN_Reference_DB/
COPY bLactam_MIC_Rscripts/ /predictor/bLactam_MIC_Rscripts/
COPY ExtractGene.pl PBP-Gene_Typer.pl pw_wrapper.sh to_json.pl transeq.pl entrypoint.sh /predictor/

RUN chmod 0755 /predictor/*.pl /predictor/*.sh

ENTRYPOINT ["/predictor/entrypoint.sh"]
