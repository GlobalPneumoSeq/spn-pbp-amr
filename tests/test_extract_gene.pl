#!/usr/bin/env perl
use strict;
use warnings;
use Test::More;

require "./ExtractGene.pl";

my $records = { contig => 'ACGT' };

is(
    extractAssemblyInterval($records, 'contig', 0, 4, 0),
    ">contig:0-4\nACGT\n",
    'forward interval retains the bedtools FASTA representation',
);
is(
    extractAssemblyInterval($records, 'contig', 0, 4, 1),
    ">contig:0-4\nACGT",
    'reverse interval retains the historical no-final-newline representation',
);
is(
    extractAssemblyInterval($records, 'contig', 0, 5, 0),
    '',
    'an interval extending beyond the contig is skipped rather than truncated',
);
is(
    extractAssemblyInterval($records, 'contig', -1, 3, 0),
    '',
    'an interval beginning before the contig is skipped',
);

done_testing();
