# PBP reference data

The three `SPN_bLactam_*-DB.faa` PBP allele FASTAs were sourced from
[`Spn_Scripts_Reference` commit `39ebc288`](https://github.com/BenJamesMetcalf/Spn_Scripts_Reference/commit/39ebc2888f1c6166ca95022050bd53b4e2435768),
the 2026-08-31 PBP database update.

`Blast_bLactam_*_prot_DB*` files are deliberately excluded: the Docker image
generates them from these FASTAs using its pinned BLAST+ release. This prevents
the source allele data and BLAST indexes from becoming inconsistent.
