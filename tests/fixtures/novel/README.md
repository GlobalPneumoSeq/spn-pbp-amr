# Novel-PBP regression fixture

`run_regression.sh` derives this fixture from `SAMN20398278.fasta` by changing
base zero-based position `200991` of contig
`NODE_1_length_271402_cov_32.208704` to `A`. This substitutes PBP2B's first
amino acid, producing the expected `2:NEW:0` allele status while retaining a
valid predictor result. Keeping the large source assembly deduplicated avoids
committing a second near-identical 2 MB fixture.
