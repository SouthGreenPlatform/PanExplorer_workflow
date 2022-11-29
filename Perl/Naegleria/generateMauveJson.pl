#!/usr/bin/perl

use strict;

my $chrom_focus = $ARGV[0];

for (my $i=37; $i >= 1; $i--){print "NL$i;";}exit;

print "[\n";
my %gene_positions;
open(G,"Map_annotation_ID_NFGwada.gff3.gff3");
while(<G>){
	my $line = $_;
	$line =~s/\n//g;$line =~s/\r//g;
	my @infos = split("\t",$line);
	if ($infos[2] eq "mRNA" && /ID=([^;]+);/){
		my $gene = $1;
		my $start = $infos[3];
		my $end = $infos[4];
		my $chr = $infos[0];
		#if ($chrom_focus && $chrom_focus ne $chr){next;}
		$gene_positions{$gene} = "$chr-$start-$end";
	}
}
close(G);

open(G,"Map_annotation_ID_NLova7.gff3.gff3");
while(<G>){
        my $line = $_;
        $line =~s/\n//g;$line =~s/\r//g;
        my @infos = split("\t",$line);
        if ($infos[2] eq "mRNA" && /ID=([^;]+);/){
                my $gene = $1;
                my $start = $infos[3];
                my $end = $infos[4];
                my $chr = $infos[0];
		#if ($chrom_focus && $chrom_focus ne $chr){next;}
                $gene_positions{$gene} = "$chr-$start-$end";
        }
}
close(G);

my $lines = "";
open(F,"orthofinder_matrix.txt");
<F>;
while(<F>){
	my $line = $_;
	$line =~s/\n//g;$line =~s/\r//g;
	my @infos = split("\t",$line);
	my $nb_found = 0;
	my $index = 0;
	for (my $i = 1; $i <= $#infos; $i++){
		my $val = $infos[$i];
		if ($val =~/\w+/){
			$nb_found++;
			$index = $i;
		}
	}
	if ($nb_found == 1){
		#print "$index\n";
	}
	#next;
	if ($nb_found == $#infos){
		my $gene1 = $infos[1];
		my $gene2 = $infos[7];
		if ($gene1 !~/,/ && $gene2 !~/,/){
			my ($chr1,$start1,$end1) = split(/-/,$gene_positions{$gene1});
			my ($chr2,$start2,$end2) = split(/-/,$gene_positions{$gene2});
			$chr2 = lc($chr2);
			#if ($chr1 eq $chr2 && $chr1 eq "$chrom_focus"){
			if ($chr1 ne $chr2 && $chr2 =~/\w+/){
			#if ($chr1 eq $chr2){
				if ($chr1=~/chr(\d+)/){
					#my $nb = $1 * 1400000;
					my $nb = 0;
					$start1 += $nb;
					$end1 += $nb;
				}
				if ($chr2=~/chr(\d+)/){
                                        #my $nb = $1 * 1400000;
					my $nb = 0;
                                        $start2 += $nb;
                                        $end2 += $nb;
                                }	
			#[{"name": "574556.4.fasta","start": 717247,"end": 718620,"strand": "-","lcb_idx": 1},
				#$lines .= "[{\"name\": \"1.fasta\",\"start\": $start1,\"end\": $end1,\"strand\": \"-\",\"lcb_idx\": 1},{\"name\": \"2.fasta\",\"start\": $start2,\"end\": $end2,\"strand\": \"-\",\"lcb_idx\": 2}],\n";
				$chr1 =~s/chr/NF/g;
				$chr2 =~s/chr/NL/g;
				print "$chr1 $start1 $end1 $chr2 $start2 $end2\n";
			}
		}
	}
}
close(F);

chop($lines);
chop($lines);
#print "$lines\n]";
