# Implementation plan: standalone pneumococcal beta-lactam MIC predictor

## Summary

Modernise and maintain a small Dockerised predictor for *Streptococcus
pneumoniae* beta-lactam MIC and SIR inference from an assembled-genome FASTA.
The initial release will recover the PBP1A, PBP2B, and PBP2X regions, generate
the existing Random Forest predictions, and write the established JSON result.

The work is deliberately a packaging and maintainability exercise, not a model
or reference-data update. It will remove avoidable dependencies (Bedtools,
EMBOSS, and Elastic Net/glmnet), replace cluster-era assumptions with pinned
container dependencies, and preserve validated Random Forest output.

## The two code bases

`spn-resistance-pbp` is the implementation target. It is an earlier,
assembly-FASTA-focused extraction of the beta-lactam component, already shaped
as a Docker application with a small stdin-FASTA-to-JSON interface. It includes
the assembly PBP extractor (`ExtractGene.pl`), the PBP allele typer
(`PBP-Gene_Typer.pl`), the MIC scripts, and an existing Dockerfile. Its
container and dependencies are dated, but its narrower scope makes it the right
place for this work.

`spn-resistance-pbp2` is the newer, broader source repository used as a
behavioural reference. It includes the original cluster-oriented full typing
pipeline and more recent changes to the MIC scripts. In particular, it contains
the R-version-safe BLOSUM62 handling for novel amino acids and reflects the
current Random Forest-only behaviour. It is not the implementation target:
its read-processing, serotyping, MLST, and grid-execution code are outside this
project's scope.

Both repositories contain versions of the MIC scripts and model/reference
files. Database and model updates are explicitly excluded from this plan. Code
changes should be compared against `spn-resistance-pbp2` where it represents a
newer intended behaviour, then validated with fixed regression fixtures in
`spn-resistance-pbp`.

## Scope

Produce a modern, reproducible Docker image for assembly FASTA input only. It
will extract PBP1A, PBP2B, and PBP2X; generate Random Forest MIC/SIR
predictions; and emit the existing JSON result contract.

The model and reference data are out of scope: retain the current supplied
database/model snapshot unchanged.

## 1. Establish a regression baseline

Create a small versioned test fixture set containing assemblies with:

- Complete forward- and reverse-strand PBP regions.
- Partial or missing PBP regions.
- A divergent/novel PBP amino-acid residue.
- At least one known successful real-world assembly.

For each fixture, record:

- Extracted PBP nucleotide sequences.
- Translated PBP amino-acid sequences.
- PBP allele status (`known`, `NEW`, `NF`, or `ERROR`).
- Final Random Forest MIC/SIR JSON.

Goal: all later changes have an explicit expected-output oracle.

## 2. Remove Bedtools; retain and update BLAST+

Replace the two `bedtools getfasta` calls in `ExtractGene.pl` with native Perl
logic that:

- Reads the matched assembly contig by ID.
- Extracts the same zero-based, half-open interval currently written to BED.
- Reverse-complements reverse-strand extracts.
- Preserves the current BLAST-derived coordinate calculations, hit ranking, and
  thresholds.

Keep `blastn` for PBP-region identification and `blastp` for existing PBP
allele identification. Update BLAST+ to a maintained, pinned release.

Goal: all baseline extracted nucleotide sequences and final RF MIC/SIR results
are identical before and after removal of Bedtools. The image contains no
`bedtools` executable or package.

## 3. Keep FASTQ support as documented future work

Do not implement FASTQ input.

Document the historical
[`PBP-Gene_Typer.pl`](https://github.com/BenJamesMetcalf/JanOw_Dependencies/blob/master/PBP-Gene_Typer.pl)
as a possible future path: it uses paired reads with `LoTrac_target.pl` to
recover PBP targets, translates them, then assigns known/new allele status with
BLASTP.

Goal: the README clearly states that the released image accepts assembly FASTA
only, and identifies the historical FASTQ workflow's extra dependency and
validation requirements.

## 4. Remove Elastic Net prediction

Delete Elastic Net prediction and classification code from the R predictor.

- Remove `glmnet` from R dependency installation.
- Do not load Elastic Net model files.
- Remove Elastic Net columns from intermediate prediction output.
- Retain only Random Forest MIC prediction and downstream MIC/SIR formatting.

Goal: source code and built image have no `glmnet` dependency or Elastic Net
model loads; baseline final MIC/SIR JSON remains unchanged.

## 5. Update novel-amino-acid substitution

Retain the newer numeric BLOSUM62 ranking logic from `spn-resistance-pbp2`.

- Convert BLOSUM scores to numeric before ranking.
- Choose the highest-scoring valid training-level amino acid.
- Make ties deterministic by retaining the first model-level amino acid with the
  highest score.

Goal: automated tests cover known residues, novel residues, ambiguous residues,
gaps, and a tie; results are stable under the pinned R version.

## 6. Modernise the container

Replace the Ubuntu 18.04 image and legacy runtime downloads with a supported
base image and pinned dependencies.

Required runtime dependencies:

- Perl and JSON support.
- Current pinned BLAST+.
- clustal-omega.
- R plus `Biostrings`, `randomForest`, `iterators`, and `foreach`.
- Bash and standard POSIX utilities.

Excluded dependencies:

- Bedtools.
- EMBOSS.
- `glmnet`.
- Grid modules, Singularity, and all `/scicomp` paths.

Goal: `docker build` succeeds from a clean checkout; `docker run --rm -i IMAGE
< assembly.fa` runs without downloads or external services.

## 7. Release verification and documentation

Add automated checks for:

- Docker build.
- Valid FASTA produces valid JSON.
- Invalid/empty FASTA fails clearly and non-zero.
- All regression fixtures match expected PBP extraction and final JSON.
- `docker image inspect` confirms the fixed image version label.
- README documents build, stdin FASTA input, JSON output, limitations, and
  reproducibility versions.

Definition of done: the image builds reproducibly, has no
Bedtools/EMBOSS/Elastic Net dependency, preserves baseline Random Forest
output, and has a passing fixture-based regression suite.
