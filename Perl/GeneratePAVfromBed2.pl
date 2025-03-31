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
                        my $genelength = $end-$start;
                        $current_gene_length = $genelength;
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
                        my $start = $infos[3];
                        my $end = $infos[4];
                        my $genelength = $end-$start;
                        $gene_lengths{$id}{$gene} += $genelength;
                        if ($line =~/protein_id=([^;]+);/){
                                $gene = $1;
                                $gene_lengths{$id}{$gene} += $genelength;
                        }
                }
        }
        close(GFF);
}
my $datestring = localtime();
print "Local date and time $datestring\n";

my %hash;
my %genes;
my $graph = Graph::Undirected->new;
my $num = 0;
foreach my $id(keys(%strains)){

        my $file = `ls $bed_directory/$id.*.bed`;
        chomp($file);
	#print "$file\n";
        open(B,$file);
        while(my $line = <B>){
                chomp($line);
                my @infos = split(/\t/,$line);
		my $node = $infos[0];
		my $start = $infos[1];
		my $end = $infos[2];
                my $gene = $infos[3];
		for (my $i = $start; $i <= $end; $i++){
			$hash{"$node:$i"} .= "$id!$gene,";
		}
        }
        close(B);
}
print "ok";
my $datestring = localtime();
print "Local date and time $datestring\n";

my %overlaps;
open(OUT2,">$PAV.position_in_graph.txt");
foreach my $position_in_graph(keys(%hash)){
	my $concatenated_genes = $hash{$position_in_graph};
	my @genes = split(/,/,$concatenated_genes);
	if (scalar @genes > 1){
		for (my $k = 1; $k <= $#genes; $k++){
			my $gene = $genes[$k];
			my $previous_gene = $genes[$k-1];
			print OUT2 "$position_in_graph $gene $previous_gene\n";
			$overlaps{$position_in_graph}{$gene} = 1;
		       	$overlaps{$position_in_graph}{$previous_gene} = 1;	
			#$graph->add_edge($gene,$previous_gene);
		}

	}
}

#foreach my $position_in_graph(keys(%overlaps)){
#	my $refhash = $overlaps{$position_in_graph};
#	my %subhash = %$refhash;
#	my @genes = keys(%subhash);
#	print OUT2 "$position_in_graph ".join(",",@genes)."\n";
#	foreach my $gene(@genes){
#		
#	}
#}
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
