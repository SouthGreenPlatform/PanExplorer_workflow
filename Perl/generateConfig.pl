#!/usr/bin/perl

use strict;
use warnings;
use YAML qw(DumpFile);

use Cwd;
my $dir = getcwd;


my $zip = $ARGV[0];
my $list = $ARGV[1];
my $out = $ARGV[2];
my $zip_fasta = $ARGV[3];

$list =~s/NZ_//g;
my @list_ids = split(/,/,$list);

my %data = ();
# case list of accessions
if ($list ne 'None'){
	foreach my $id(@list_ids){
		push @{$data{"ids"}}, "$id";
	}
}
# case fasta+gff
elsif ($zip ne "None" && $zip_fasta ne "None"){
	system("rm -rf $zip.genomeszip");
	mkdir("$zip.genomeszip");
	chdir("$zip.genomeszip");
	system("cp -rf $zip ./genomes.zip");
	system("cp -rf $zip_fasta ./fasta_genomes.zip");
	system("unzip fasta_genomes.zip");
	system("unzip genomes.zip");
	unlink("fasta_genomes.zip");
	unlink("genomes.zip");
	open(LS,"ls *gff |");
	while(my $line = <LS>){
		chomp($line);
		my $concat = "$zip.genomeszip/".$line;
		if ($line =~/(.*)\.gff/){
			my $name = $1;
			if ($concat =~/\w+/){
				$data{"input_genomes"}{"$name"}{"gff3"} = "$concat";
				$data{"input_genomes"}{"$name"}{"name"} = "$name";
			}
		}
	}
	close(LS);

	open(LS,"ls *fasta |");
        while(my $line = <LS>){
                chomp($line);
                my $concat = "$zip.genomeszip/".$line;
                if ($line =~/(.*)\.fasta/){
                        my $name = $1;
                        if ($concat =~/\w+/){
				$data{"input_genomes"}{"$name"}{"fasta"} = "$concat";
                        }
                }
        }
        close(LS);

}
elsif ($zip ne "None"){
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
