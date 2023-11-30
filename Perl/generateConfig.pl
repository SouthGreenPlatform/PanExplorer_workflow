#!/usr/bin/perl

use strict;
use warnings;
use YAML qw(DumpFile);

use Cwd;
my $dir = getcwd;


my $zip = $ARGV[0];
my $list = $ARGV[1];
my $out = $ARGV[2];

$list =~s/NZ_//g;
my @list_ids = split(/,/,$list);

my %data = ();
foreach my $id(@list_ids){
	push @{$data{"ids"}}, "$id";
}
if ($zip ne "None"){
	system("rm -rf $zip.genomeszip");
	mkdir("$zip.genomeszip");
	chdir("$zip.genomeszip");
	my $head = `head -1 $zip`;
	if ($head =~/LOCUS/){
		push @{$data{"input_genbanks"}}, "$zip";
	}
	else{
		system("cp -rf $zip ./genomes.zip");
		system("unzip genomes.zip");
		unlink("genomes.zip");
		open(LS,"ls |");
		while(my $line = <LS>){
			chomp($line);
			my $concat = "$zip.genomeszip/".$line;
			if ($concat =~/\w+/){
				push @{$data{"input_genbanks"}}, "$concat";
			}
		}
		close(LS);
	}
}
#else{
#	push @{$data{"input_genbanks"}}, "";
#}
chdir($dir);

DumpFile($out, \%data);
