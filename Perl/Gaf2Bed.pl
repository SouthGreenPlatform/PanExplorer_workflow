#!/usr/bin/perl

use strict;


my $gaf = $ARGV[0];
my $outbed = $ARGV[1];

open(O,">$outbed");

my %segments_chr;
my %segments_of_gene;
my %segments;
my %links;
my $max = 0;
open(G,$gaf);
while(<G>){
	my $line = $_;
	$line =~s/\n//g;$line =~s/\r//g;
	my @infos = split(/\t/,$line);
	my $gene = $infos[0];
	my $path = $infos[5];
	my $start_in_path = $infos[7];
	my $end_in_path = $infos[8];
	print O "$path	$start_in_path	$end_in_path	$gene\n";
}
close(G);
close(O);
