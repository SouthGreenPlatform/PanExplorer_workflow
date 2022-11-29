#!/usr/bin/perl

use strict;

my $gff = $ARGV[0];
my $window_size = $ARGV[1];

my %counts;
my %hash;
open(F,$gff);
while(<F>){
	if (/^#/){next;}
	my @infos = split("\t",$_);
	my $chr = $infos[0];
	#$chr =~s/chr//g;
	my $source = $infos[1];
	my $feature = $infos[2];
	my $start = $infos[3];
	my $end = $infos[4];
	my $window_num = int($start/$window_size)+1;
	my $window_num2 = int($end/$window_size)+1;
	if ($source eq "repeatmasker" && $feature eq "match"){
		$feature = "repeats";
	}
	if ($feature ne "gene" && $feature ne "CDS" && $feature ne "repeats" && $feature !~/UTR/ && $feature ne "exon"){next;}
	$counts{$chr}{$window_num}{$feature}++;
	if ($window_num == $window_num2){
		for (my $i = $start; $i <= $end; $i++){
			$hash{$chr}{$window_num}{$feature}{$i} = 1;
		}
	}
	else{
		for (my $i = $start; $i <= ($window_num*$window_size); $i++){
			$hash{$chr}{$window_num}{$feature}{$i} = 1;
                }
		for (my $j = ($window_num2-1)*$window_size; $j <= $end; $j++){
			$hash{$chr}{$window_num2}{$feature}{$j} = 1;
                }
		if (($window_num2 - $window_num) > 1){
			for (my $window_num_new = $window_num + 1; $window_num_new <= $window_num2 - 1; $window_num_new++){
				for (my $j = (($window_num_new-1)*$window_size); $j <= ($window_num_new*$window_size); $j++){
					$hash{$chr}{$window_num_new}{$feature}{$j} = 1;
				}
			}
		}
		
	}
}
close(F);

#my $refhash2 = $hash{"Chr13"}{"1"}{"exon"};
#my %subhash2 = %$refhash2;
#print scalar %subhash2;
#exit;

open(G2,">gene_counts.txt");
open(G,">gene_density.txt");
open(I,">exon_density.txt");
open(R,">repeat_density.txt");
foreach my $chr(sort {$a<=>$b} keys(%hash)){
	#print G "$chr	0	1	0\n";
	#print G "$chr	1	2	1\n";
	#print I "$chr	0	1	0\n";
	#print I "$chr	0	1	0\n";
	my $refhash = $hash{$chr};
	my %subhash = %$refhash;
	foreach my $window_num(sort {$a<=>$b} keys(%subhash)){

		# genes
		my $proportion_gene = 0;
		if ($hash{$chr}{$window_num}{"gene"}){
			my $refhash2 = $hash{$chr}{$window_num}{"gene"};
			my %subhash2 = %$refhash2;
			$proportion_gene = scalar %subhash2 / $window_size;
			if ($proportion_gene > 1.1){
				#print "$chr $window_num $proportion_gene ".scalar %subhash2."\n";exit;
			}
		}

		# exon
		my $proportion_exon = 0;
		if ($hash{$chr}{$window_num}{"CDS"}){
			my $refhash2 = $hash{$chr}{$window_num}{"CDS"};
			my %subhash2 = %$refhash2;
			$proportion_exon = scalar %subhash2 / $window_size;
		}

		# repeat
		my $proportion_repeat = 0;
                if ($hash{$chr}{$window_num}{"repeats"}){
                        my $refhash2 = $hash{$chr}{$window_num}{"repeats"};
                        my %subhash2 = %$refhash2;
                        $proportion_repeat = scalar %subhash2 / $window_size;
			#$proportion_repeat = $proportion_repeat + $proportion_gene;
                }
		else{
			#$proportion_repeat = $proportion_gene;
			$proportion_repeat = 0;
		}
		my $start = ($window_num-1)*$window_size;
		my $end = $window_num*$window_size;
		print I "$chr	$start	$end	$proportion_exon\n";
		print G "$chr	$start	$end	$proportion_gene\n";
		print R "$chr	$start	$end	$proportion_repeat\n";
		if ($counts{$chr}{$window_num}{"gene"}){
			print G2 "$chr	$start	$end	".$counts{$chr}{$window_num}{"gene"}."\n";
		}
		else{
			print G2 "$chr	$start	$end	0\n";
		}
	}
}
close(G);
close(I);
close(R);
close(G2);

