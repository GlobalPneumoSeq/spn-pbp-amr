#!/usr/bin/env perl
use strict;
use warnings;
use JSON;

my %data;
print STDERR "Parsing output to JSON.";

my $line = <>;
die "Predictor output did not contain a result row\n"
    unless defined $line && $line =~ /\S/;
die "Predictor output contained more than one result row\n" if <>;
my @fields = split /\s+/, $line;

die "Predictor output did not contain a complete result row\n"
    unless @fields >= 25;

$data{"pbp1a"} = $fields[1];
$data{"pbp2b"} = $fields[2];
$data{"pbp2x"} = $fields[3];
$data{"penMic"} = $fields[4].' '.$fields[5];
$data{"penMeningitis"} = $fields[6];
$data{"penNonMeningitis"} = $fields[7];
$data{"amxMic"} = $fields[8].' '.$fields[9];
$data{"amx"} = $fields[10];
$data{"memMic"} = $fields[11].' '.$fields[12];
$data{"mem"} = $fields[13];
$data{"ctxMic"} = $fields[14].' '.$fields[15];
$data{"ctxMeningitis"} = $fields[16];
$data{"ctxNonMeningitis"} = $fields[17];
$data{"croMic"} = $fields[18].' '.$fields[19];
$data{"croMeningitis"} = $fields[20];
$data{"croNonMeningitis"} = $fields[21];
$data{"cxmMic"} = $fields[22].' '.$fields[23];
$data{"cxm"} = $fields[24];

my $str = encode_json \%data;
print STDOUT $str;
