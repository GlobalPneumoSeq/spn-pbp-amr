# Error-PBP regression fixture

`run_regression.sh` derives this fixture from `examples/SAMN20398278.fasta` by
inserting `A` at zero-based position `200995` of contig
`NODE_1_length_271402_cov_32.208704`. The resulting PBP2B frameshift remains
detectable at the nucleotide level but cannot be assigned an allele, producing
the expected `2:ERROR:0` status and a valid JSON predictor result.
