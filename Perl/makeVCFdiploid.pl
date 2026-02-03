#!/usr/bin/perl

use strict;

my $in = $ARGV[0];
my $out = $ARGV[1];
my $ref = $ARGV[2];

open(I,$in);
open(O,">$out");
while(<I>){
	if (!/^#/){
		my $line = $_;
		$line =~s/\n//g;$line =~s/\r//g;
		my @infos = split(/\t/,$line);
		print O "ref";
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
