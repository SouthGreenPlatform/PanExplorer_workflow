#!/usr/bin/perl

use strict;

use Graph::Undirected;

my $strain_info = $ARGV[0];
my $bed_directory = $ARGV[1];
my $PAV = $ARGV[2];
my $threshold_coverage = 30;

my %strains;
open(S,$strain_info);
while (my $line = <S>) {
	chomp $line;
	my ($id,$name) = split(/\t/,$line);
	$strains{$id} = $name;
}
close(S);

my %gene_lengths;
foreach my $id(keys(%strains)){
	my $f = "$bed_directory/$id.gene_segments.gene_length.txt";
	open(F,$f);
	while(my $line =<F>){
		chomp($line);
		my ($gene,$length) = split(/\t/,$line);
		$gene_lengths{$id}{$gene} = $length;
	}
	close(F);	
}


my %genes;
my $graph = Graph::Undirected->new;
my $num = 0;
foreach my $id(keys(%strains)){
	my $file = "$bed_directory/$id.gene_segments.bed";
	open(B,$file);
	while(<B>){
		my $line = $_;
		$line =~s/\n//g;$line =~s/\r//g;
		my @infos = split(/\t/,$line);
		my $gene1 = $infos[3];
		$genes{"$id:$gene1"} = 1;
	}
	close(B);

	foreach my $id2(keys(%strains)){
		my $file2 = "$bed_directory/$id2.gene_segments.bed";
		if (-e $file && -e $file2 && $file ne $file2){
			$num++;
			system("bedtools intersect -a $file -b $file2 -wo >$PAV.$num.intersect.out");
			my %cumul_match;
			my %cumul_match2;
			open(INTER,"$PAV.$num.intersect.out");
			while(<INTER>){
				my $line = $_;
				$line =~s/\n//g;$line =~s/\r//g;
				my ($segment,$start_gene1,$end_gene1,$gene1,$segment,$start_gene2,$end_gene2,$gene2,$size_overlap) = split(/\t/,$line);
				$cumul_match{"$gene1;$gene2"}+=$size_overlap;
			}
			close(INTER);

			foreach my $pair(keys(%cumul_match)){
				my $size1 = $cumul_match{$pair};
				my ($gene1,$gene2) = split(/;/,$pair);
				my $gene1_length = $gene_lengths{$id}{$gene1};
				my $gene2_length = $gene_lengths{$id2}{$gene2};
				my $size_match_gene = $cumul_match{"$gene1;$gene2"};
				my $percentage_overlap_gene1 = ($size_match_gene/$gene1_length)*100;
				my $percentage_overlap_gene2 = ($size_match_gene/$gene2_length)*100;
				#print "$pair $size_match_gene $percentage_overlap_gene1 $percentage_overlap_gene2\n";
				if ($percentage_overlap_gene1 > $threshold_coverage && $percentage_overlap_gene2 > $threshold_coverage){
					$graph->add_edge("$id:$gene1","$id2:$gene2");
					delete($genes{"$id:$gene1"});
					delete($genes{"$id2:$gene2"});
				}
			}
		}
	}
}

open(OUT,">$PAV");
print OUT "ClutserID";
foreach my $id(sort keys(%strains)){
	my $name = $strains{$id};
	print OUT "\t".$name;
}
print OUT "\n";
my $clnum = 0;
my @cc = $graph->connected_components();
foreach my $component (@cc){
	$clnum++;
	print OUT $clnum;
        my @genes = @$component;
	my %h;
	foreach my $gene(@genes){
		my ($id,$genename) = split(/:/,$gene);
		$h{$id}.="$genename,";
	}	
	foreach my $id(sort keys(%strains)){
		my $ids = "-";
		if ($h{$id}){
			$ids = $h{$id};	
			chop($ids);
		}
		print OUT "\t".$ids;
	}
	print OUT "\n";
}
#######################
# add singletons
#######################
foreach my $gene(keys(%genes)){
	$clnum++;
	print OUT $clnum;
	my ($id,$genename) = split(/:/,$gene);
	my %h;
	$h{$id}.="$genename,";
	foreach my $id(sort keys(%strains)){
                my $ids = "-";
                if ($h{$id}){
                        $ids = $h{$id};
                        chop($ids);
                }
                print OUT "\t".$ids;
        }
        print OUT "\n";	
}
close(OUT);
