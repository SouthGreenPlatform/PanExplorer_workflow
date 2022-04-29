#!/usr/bin/perl

use strict;
use File::Basename;
my $dirname = dirname(__FILE__);

my $pav_matrix = $ARGV[0];
my $prot_dir = $ARGV[1];
my $cog_outfile = $ARGV[2];
my $cog_stats = $ARGV[3];
my $cog_stats2 = $ARGV[4];
my $strain_info_file = $ARGV[5];


my %cogs_categories = (
        "J"=>"INFORMATION STORAGE AND PROCESSING",
        "A"=>"INFORMATION STORAGE AND PROCESSING",
        "K"=>"INFORMATION STORAGE AND PROCESSING",
        "L"=>"INFORMATION STORAGE AND PROCESSING",
        "B"=>"INFORMATION STORAGE AND PROCESSING",
        "D"=>"CELLULAR PROCESSES AND SIGNALING",
        "Y"=>"CELLULAR PROCESSES AND SIGNALING",
        "V"=>"CELLULAR PROCESSES AND SIGNALING",
        "T"=>"CELLULAR PROCESSES AND SIGNALING",
        "M"=>"CELLULAR PROCESSES AND SIGNALING",
        "N"=>"CELLULAR PROCESSES AND SIGNALING",
        "Z"=>"CELLULAR PROCESSES AND SIGNALING",
        "W"=>"CELLULAR PROCESSES AND SIGNALING",
        "U"=>"CELLULAR PROCESSES AND SIGNALING",
        "O"=>"CELLULAR PROCESSES AND SIGNALING",
        "C"=>"METABOLISM",
        "G"=>"METABOLISM",
        "E"=>"METABOLISM",
        "F"=>"METABOLISM",
        "H"=>"METABOLISM",
        "I"=>"METABOLISM",
        "P"=>"METABOLISM",
        "Q"=>"METABOLISM",
        "R"=>"POORLY CHARACTERIZED",
        "S"=>"POORLY CHARACTERIZED"
);

my %strain_names;
open(S,$strain_info_file);
while(<S>){
	my $line = $_;
	$line =~s/\n//g;$line =~s/\r//g;
	my ($id,$strain_name) = split(/\t/,$line);
	$strain_names{$id} = $strain_name;
}
close(S);

my %strain_of_prot;
my %proteins;
open( DIR, "ls $prot_dir/*pep |" );
while(<DIR>) {
	my $filename = $_;
	my $strain;
	my $id;
	if ($filename =~/\/([^\/]*).pep/){
		$strain = $1;
	}
	#open(F,"zcat $filename|" );
	open(F,"$filename" );
	while(<F>){
		if (/>(.*)/){
			$id = $1;
			$strain_of_prot{$id} = $strain;
		}
		else{
			$proteins{$id}.= $_;
		}
	}
	close(F);
}
closedir(DIR);


my %genes_of_cluster;
open(S,">$pav_matrix.selection_prot.fa");
open(O,$pav_matrix);
<O>;
while(<O>){
        my $line = $_;
        $line =~s/\n//g;$line =~s/\r//g;
        my @infos = split(/\t/,$line);
        my $cluster = $infos[0];
        print S ">$cluster\n";
	for (my $i=1; $i <= $#infos; $i++){
		my @genenames = split(/,/,$infos[$i]);
		foreach my $genename(@genenames){
                	if ($genename=~/\w+/){
                        	$genes_of_cluster{$cluster} .= "$genename,";
			}
		}
	}
        for (my $i=1; $i <= $#infos; $i++){
		my @genenames = split(/,/,$infos[$i]);
		foreach my $genename(@genenames){
			if ($genename=~/\w+/){
				print S $proteins{$genename};
				last;
			}
		}
        }
}
close(O);
close(S);

system("rpsblast+ -query $pav_matrix.selection_prot.fa -db $dirname/../COG/Cog -out $pav_matrix.selection.rps-blast.out -evalue 1e-2 -outfmt '7 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs'");

system("perl $dirname/../COG/bac-genomics-scripts/cdd2cog/cdd2cog.pl -r $pav_matrix.selection.rps-blast.out -c $dirname/../COG/cddid.tbl -f $dirname/../COG/fun.txt -w $dirname/../COG/whog -a");

open(COG,">$cog_outfile");
my %cogs;
my %count_letter;
open(C,"results/protein-id_cog.txt");
while(<C>){
	my $line = $_;
	$line =~s/\n//g;$line =~s/\r//g;
	my @infos = split(/\t/,$line);
	my $coginfo = "";
	my $cluster = $infos[0];
	for (my $i = 1; $i <= $#infos; $i++){
		$coginfo .= "\t".$infos[$i];
	}
	my $gene_list = $genes_of_cluster{$cluster};
	chop($gene_list);
	my @genes = split(",",$gene_list);
	foreach my $gene(@genes){
		my @letters = split(/\t/,$coginfo);
		my $strain = $strain_of_prot{$gene};
		for (my $i = 2; $i <= $#letters; $i++){
			my $letter = $letters[$i];
			my $cat = $cogs_categories{$letter};
			$count_letter{$strain}{$letter}++;
			$count_letter{$strain}{$cat}++;
		}
		print COG $gene.$coginfo."\n";
	}
}
close(C);
system("rm -rf results");
close(COG);

my @cat_of_cat = ("INFORMATION STORAGE AND PROCESSING","METABOLISM","CELLULAR PROCESSES AND SIGNALING","POORLY CHARACTERIZED");
my @cog_categories = ("D","M","N","O","T","U","V","W","Y","Z","A","B","J","K","L","C","E","F","G","H","I","P","Q","R","S");
open(COG_STAT,">$cog_stats");
open(COG_STAT2,">$cog_stats2");
print COG_STAT "COG\t".join("\t",@cog_categories)."\n";
print COG_STAT2 "COG\t".join("\t",@cat_of_cat)."\n";
foreach my $strain(keys(%count_letter)){
	my $strain_name = $strain_names{$strain};
	print COG_STAT "$strain_name";
	print COG_STAT2 "$strain_name";
	my $ref_subhash = $count_letter{$strain};
	my %subhash = %$ref_subhash;
	foreach my $letter(@cog_categories){
		my $n = 0;
		if ($count_letter{$strain}{$letter}){$n = $count_letter{$strain}{$letter};}
		print COG_STAT "\t".$n;
	}
	print COG_STAT "\n";

	foreach my $cat(@cat_of_cat){
		my $n = 0;
		if ($count_letter{$strain}{$cat}){$n = $count_letter{$strain}{$cat};}			
		print COG_STAT2 "\t".$n;
	}
	print COG_STAT2 "\n";
}
close(COG_STAT);
close(COG_STAT2);

