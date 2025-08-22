#!/usr/bin/perl

use strict;

my $in = $ARGV[0];
my $out = $ARGV[1];

open(IN,$in);
open(OUT,">$out");
my $id = "";
my $n = 0;
while(<IN>){
	my $line = $_;
	$line =~s/\n//g;$line =~s/\r//g;
	if ($line =~/Sequence: (.*)$/){$id=$1;}
	elsif ($line =~/^(\d+) (\d+)/){
		$n++;
		my @infos = split(/ /,$line);
		print OUT "$id	$infos[0]	$infos[1]	repeat".$n."_".$infos[13]."_".$infos[3]."\n";
	}
}
close(IN);
close(OUT);
