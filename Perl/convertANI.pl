#!/usr/bin/perl

use strict;

my $file = $ARGV[0];
my $metadata = $ARGV[1];


my %strains;
open(F,$metadata);
while(<F>){
	my $line = $_;
	$line =~s/\n//g;$line =~s/\r//g;
	my ($path_genome,$strain) = split(/\t/,$line);
	$strains{$path_genome} = $strain;
}
close(F);

my %ANIs;
my %genomes;
my %genome_names;
my $num_line = 0;
open(F,$file);
<F>;
while(<F>){
	$num_line++;	
	my $line = $_;
	$line =~s/\n//g;$line =~s/\r//g;
	my @infos = split(/\t/,$line);
	my $genome = $infos[0];
	$genome = $strains{$genome};
	$genome_names{$genome}++;
	if ($genome_names{$genome} > 1){
		$genome = $genome . ".". $genome_names{$genome};
	}

	$genomes{$num_line} = $genome;
	for (my $i = 1; $i <= $#infos; $i++){
		$ANIs{$i}{$num_line} = $infos[$i];
		$ANIs{$num_line}{$i} = $infos[$i];
	}
}
close(F);

print "Genomes";

foreach my $i(sort keys(%ANIs)){
	print "\t".$genomes{$i};
}
print "\n";
foreach my $i(sort keys(%ANIs)){
	print $genomes{$i};
	foreach my $j(sort keys(%ANIs)){
		my $ani = $ANIs{$i}{$j};
		if ($ani eq "NA"){$ani = 0;}
		if (!$ani && $i == $j){
			$ani = 100;
		}
		print "\t".$ani;
	}
	print "\n";
}
