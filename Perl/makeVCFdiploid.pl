#!/usr/bin/perl

use strict;

my $in = $ARGV[0];
my $out = $ARGV[1];
my $ref = $ARGV[2];
# optional 4th argument "keepchrom": keep the chromosome of the reference path in the CHROM column
# (without the PanSN prefix "genome#hap#") instead of "ref" (several chromosomes: eukaryotes)
my $keepchrom = (defined $ARGV[3] && $ARGV[3] eq "keepchrom");

open(I,$in);
open(O,">$out");
while(<I>){
	if (!/^#/){
		my $line = $_;
		$line =~s/\n//g;$line =~s/\r//g;
		my @infos = split(/\t/,$line);
		if ($keepchrom){
			my $chrom = $infos[0];
			$chrom =~s/^[^#]*#[^#]*#//;
			print O $chrom;
		}
		else{
			print O "ref";
		}
		for (my $i = 1; $i <= $#infos; $i++){
			my $val = $infos[$i];
			if ($i > 8){
				print O "\t"."$val/$val";
			}
			else{
				print O "\t$val";
			}
		}
		print O "\t0/0";
		print O "\n";
	}
	elsif (/^#CHROM/){
		my $line = $_;
		$line =~s/\n//g;$line =~s/\r//g;
		print O $line;
		print O "\t$ref";
		print O "\n";
	}
	else{
		print O $_;
	}
}
close(I);
close(O);
