#!/usr/bin/perl

use strict;

my $zip = $ARGV[0];
my $out = $ARGV[1];

open(O,">$out");
if (-e $zip){
	system("rm -rf $zip.genomeszip");
	mkdir("$zip.genomeszip");
	chdir("$zip.genomeszip");
	my $head = `head -1 $zip`;
	if ($head =~/LOCUS/){
		print "Dans ce cas\n";
		print O "$zip\n";
		exit;
	}
	system("cp -rf $zip ./genomes.zip");
	system("unzip genomes.zip");
	unlink("genomes.zip");
	open(LS,"ls |");
	while(<LS>){
		print O "$zip.genomeszip/".$_;
	}
	close(LS);
	
}
close(O);

