# Standalone pneumococcal beta-lactam MIC predictor

This repository packages the assembly-FASTA beta-lactam predictor used for
*Streptococcus pneumoniae*. It extracts PBP1A, PBP2B, and PBP2X, assigns PBP
alleles, predicts Random Forest MIC values, and emits the established JSON
result contract.

The supplied model and reference database are retained unchanged. This update
modernises the runtime and removes the old Bedtools, EMBOSS, and Elastic Net
requirements.

## Build

Linux is the primary target. Docker is required.

```sh
docker build --rm -t spn_pbp_amr:0.2.0 .
```

The image records its image version, model-data snapshot (`2023-01-12`), and
BLAST+ version (`2.16.0`) as OCI labels. The base image is pinned by digest in
the Dockerfile.

## Run

The image accepts one assembled genome FASTA on standard input and writes one
JSON object to standard output:

```sh
cat genome.fa | docker run --rm -i spn_pbp_amr:0.2.0 > result.json
```

Empty or malformed FASTA input fails with a non-zero exit status and a clear
message on standard error. The input must be an assembled nucleotide FASTA;
FASTQ and raw-read input are not supported by this image.

The JSON contains the existing beta-lactam fields, including PBP allele codes,
MIC values, and meningitis/non-meningitis SIR calls. It is intentionally not a
new schema.

## Historical FASTQ path

The original `PBP-Gene_Typer.pl` workflow documented in
[JanOw_Dependencies](https://github.com/BenJamesMetcalf/JanOw_Dependencies/blob/master/PBP-Gene_Typer.pl)
used paired reads and `LoTrac_target.pl` to recover PBP targets before
translation and allele assignment. That workflow requires additional read
processing, reference, and validation dependencies and is future work here.

## Regression fixtures

`test/fixtures/always-on/` contains three assemblies selected from the
Pathogenwatch Always-on collection, their API metadata, expected JSON, and
the extracted nucleotide/translated amino-acid artifacts used by regression
checks. The metadata records the current API MLST and serotype assignments at
fixture retrieval time; it does not contain an API key.

Run the complete Linux/Docker regression suite with:

```sh
tests/run_regression.sh
```

The suite builds the image, checks BLOSUM substitution behavior, verifies
invalid-input failures, compares extracted PBP and translated sequences, and
compares canonicalized JSON output.

## Reproducibility and limitations

- The container is pinned to the `r-base:4.4.2` image digest used by the
  Dockerfile.
- R uses Bioconductor 3.20, `Biostrings`, `randomForest`, `iterators`, and
  `foreach`; Elastic Net/glmnet is not installed or loaded.
- BLAST+ 2.16.0 is verified against NCBI's published archive checksum, and
  Debian's `clustalo` package is installed at build time; the image performs
  no runtime downloads.
- The current model/reference snapshot is not retrained or updated by this
  repository.
- Predictions are research software outputs. Validate suitability for the
  intended surveillance or clinical context before use.

Please credit the original authors of the
[Spn_Scripts_Reference](https://github.com/BenJamesMetcalf/Spn_Scripts_Reference)
software in resulting publications.
