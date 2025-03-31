#!/usr/bin/perl

use strict;

use Graph::Undirected;

my $strain_info = $ARGV[0];
my $bed_directory = $ARGV[1];
my $PAV = $ARGV[2];
my $threshold_coverage = $ARGV[3];

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

        my $gff = `ls $bed_directory/$id*.gff`;
        open(GFF,$gff);
        my $current_gene_length;
        while(my $line = <GFF>){
                chomp($line);
                my @infos = split(/\t/,$line);
                if ($infos[2] eq "gene" && $line =~/ID=([^;]+);*/){
                        my $chr = $infos[0];
                        my $start = $infos[3];
                        my $end = $infos[4];
                        my $gene = $1;
			$gene =~s/gene://g;
                        my $genelength = $end-$start;
                        $current_gene_length = $genelength;
			#print "$id $gene $genelength\n";
                        $gene_lengths{$id}{$gene} = $genelength;
                }
                if ($infos[2] eq "CDS"){
                        my $gene;
                        if ($line =~/Parent=([^;]+);*/){
                                $gene = $1;
                        }
                        elsif ($line =~/ID=([^;]+);*/){
                                $gene = $1;
                        }
			$gene =~s/gene://g;
                        my $start = $infos[3];
                        my $end = $infos[4];
                        my $genelength = $end-$start;
                        $gene_lengths{$id}{$gene} += $genelength;
                        if ($line =~/protein_id=([^;]+);/){
                                $gene = $1;
				#print "$id $gene $genelength\n";
                                $gene_lengths{$id}{$gene} += $genelength;
                        }
                }
        }
        close(GFF);
}
my $datestring = localtime();
print "Local date and time $datestring\n";

my %genes_of_node;
my %overlaps;
my %genes;
my $graph = Graph::Undirected->new;
my $num = 0;
foreach my $id(sort keys(%strains)){
	$num++;
        my $file = `ls $bed_directory/$id.*.bed`;
        chomp($file);
	my $sorted_file = $file;
	$sorted_file =~s/bed/sorted.bed/g;
	system("sort -k1,1 -k2,2n $file >$sorted_file");
	open(B,$file);
        while(my $line = <B>){
                chomp($line);
                my @infos = split(/\t/,$line);
                my $node = $infos[0];
                my $start = $infos[1];
                my $end = $infos[2];
                my $gene = $infos[3];
                $genes{"$id!$gene"} = 1;
        }
        close(B);
}
my $datestring = localtime();
print "Local date and time $datestring\n";

system("bedtools multiinter -i $bed_directory/*sorted.bed >$PAV.intersection.txt");

my $datestring = localtime();
print "Local date and time $datestring\n";

my %hash;
my $num = 0;
my %asso_num;
open(LS,"ls $bed_directory/*sorted.bed |");
while(my $file = <LS>){
	$num++;
	chomp($file);
	system("bedtools intersect -a $PAV.intersection.txt -b $file -wo >$file.intersect.out");
	if ($file =~/([^\/]*)\.+\w+.sorted.bed/){
		my $id = $1;
		print "$id\n";
		$asso_num{$num} = $id;
	}
	unlink($file);
	open(F,"$file.intersect.out");
	while(<F>){
		my @infos = split(/\t/,$_);
		my $section = $infos[0]."?".$infos[1]."-".$infos[2];
		my $gene = $infos[$#infos-1];
		$genes_of_node{$section}{$num} = $gene;
	}
	close(F);
}

my $datestring = localtime();
print "Local date and time $datestring\n";

my %overlaps;
open(INTER,"$PAV.intersection.txt");
while(my $line = <INTER>){
	chomp($line);
	my @infos = split(/\t/,$line);
	my $node = $infos[0];
	my $start = $infos[1];
	my $end = $infos[2];
	my $size = $end - $start;
	my $list = $infos[4];
	my @nums = split(/,/,$list);
	if (scalar @nums > 1){
		for (my $k=1; $k <= $#nums; $k++){
			my $num = $nums[$k];
			my $previous_num = $nums[$k-1];
			my $id = $asso_num{$num};
			my $previous_id = $asso_num{$previous_num};
			my $gene = $id."!".$genes_of_node{"$node?$start-$end"}{$num};
			my $previous_gene = $previous_id."!".$genes_of_node{"$node?$start-$end"}{$previous_num};
			
			$overlaps{"$gene"."#"."$previous_gene"}+=$size;

			#$graph->add_edge($gene,$previous_gene);
			#delete($genes{"$gene"});
			#delete($genes{"$previous_gene"});
		}
	}
}
close(INTER);

open(OUT2,">$PAV.out2");
foreach my $pair(keys(%overlaps)){
	#print OUT2 "$pair\n";
	my $size_overlap = $overlaps{$pair};
	my ($gene1,$gene2) = split(/#/,$pair);
	my ($id1,$gene_1) = split(/!/,$gene1);
	my ($id2,$gene_2) = split(/!/,$gene2);
	my $gene_length1 = $gene_lengths{$id1}{$gene_1};
	my $gene_length2 = $gene_lengths{$id2}{$gene_2};

	#print "$id1 $gene_1 => $gene_length1\n";
	#print "$id2 $gene_2 => $gene_length2\n";
	
	my $percentage_overlap_gene1 = 2;
	my $percentage_overlap_gene2 = 2;
	if ($gene_length1 && $gene_length2){
		$percentage_overlap_gene1 = ($size_overlap/$gene_length1)*100;
		$percentage_overlap_gene2 = ($size_overlap/$gene_length2)*100;
	}
	if ($percentage_overlap_gene1 > $threshold_coverage && $percentage_overlap_gene2 > $threshold_coverage){
		print OUT2 "$gene1,$gene2 $size_overlap $gene_length1\n";
		$graph->add_edge($gene1,$gene2);
		delete($genes{"$gene1"});
		delete($genes{"$gene2"});
	}
}
close(OUT2);

my $datestring = localtime();
print "Local date and time $datestring\n";

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
                my ($id,$genename) = split(/!/,$gene);
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
foreach my $gene(keys(%genes)){
        $clnum++;
        print OUT $clnum;
        my ($id,$genename) = split(/!/,$gene);
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
my $datestring = localtime();
print "Local date and time $datestring\n";
