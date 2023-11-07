#!/usr/bin/perl

use strict;

use Graph::Undirected;

my $strain_info = $ARGV[0];
my $bed_directory = $ARGV[1];
my $PAV = $ARGV[2];

my %strains;
open(S,$strain_info);
while (my $line = <S>) {
	chomp $line;
	my ($id,$name) = split(/\t/,$line);
	$strains{$id} = $name;
}
close(S);


my %genes;
my $graph = Graph::Undirected->new;
my $num = 0;
foreach my $id(keys(%strains)){
	my $file = "$bed_directory/$id.nuc.pangenome.gaf.bed";
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
		my $file2 = "$bed_directory/$id2.nuc.pangenome.gaf.bed";
		if ($file ne $file2){
			$num++;
			system("bedtools intersect -a $file -b $file2 -f 0.2 -F 0.2 -wa -wb >$PAV.$num.intersect.out");
			open(INTER,"$PAV.$num.intersect.out");
			while(<INTER>){
				my $line = $_;
				$line =~s/\n//g;$line =~s/\r//g;
				my @infos = split(/\t/,$line);
				my $gene1 = $infos[3];
				my $gene2 = $infos[7];
				$graph->add_edge("$id:$gene1","$id2:$gene2");
				delete($genes{"$id:$gene1"});
				delete($genes{"$id2:$gene2"});
			}
			close(INTER);
		}
	}
}

open(OUT,">$PAV");
print OUT "ClutserID\t";
foreach my $id(sort keys(%strains)){
	my $name = $strains{$id};
	print OUT "\t".$name;
}
print OUT "\n";
my $clnum = 0;
my @cc = $graph->connected_components();
foreach my $component (@cc){
	$clnum++;
	print OUT $clnum."\t";
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
foreach my $gene(keys(%genes)){
	$clnum++;
	print OUT $clnum."\t";
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
