#!/usr/bin/perl

use strict;
use File::Basename;
my $dirname = dirname(__FILE__);

my $pav_matrix = $ARGV[0];
my $prot_dir = $ARGV[1];
my $cog_outfile = $ARGV[2];
my $cog_stats = $ARGV[3];
my $cog_stats2 = $ARGV[4];
my $cog_clusters = $ARGV[5];
my $strain_info_file = $ARGV[6];


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

my $use_func_files = `ls $prot_dir/*func | wc -l`;
$use_func_files =~s/\n//g;$use_func_files =~s/\r//g;

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

my %cogs_of_cluster;
my %functions;
my %accessory_clusters;
my %genes_of_cluster;
my %cluster_of_gene;
open(S,">$pav_matrix.selection_prot.fa");
open(O,$pav_matrix);
<O>;
my $n = 0;
while(<O>){
	$n++;
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
                if ($infos[$i] =~/-/){$accessory_clusters{$cluster} = 1;}
        }
        LINE: for (my $i=1; $i <= $#infos; $i++){
		my @genenames = split(/,/,$infos[$i]);
		foreach my $genename(@genenames){
			if ($genename=~/\w+/){
				#$genename =~s/\|/_/g;
				print S $proteins{$genename};
				$cluster_of_gene{$genename} = $cluster;
				if ($use_func_files > 0){
					my $function = `grep '$genename' $prot_dir/*func`;
        	                        $function =~s/\n//g;$function =~s/\r//g;
                	                my @t = split(/$genename/,$function);
                        	        if ($t[1] =~/-\s+(.*)/){
                                	        $function = $1;
                                        	$function =~s/'//g;
                                	}
	                                $functions{$cluster} = $function;
					last LINE;
				}
			}
		}
        }
}
close(O);
close(S);


#########################################################################################################################################
# Test if COG are already provided in Genbank files. If provided, rpsblast is not launched, we get the COG information from Genbank files
#########################################################################################################################################
my %cog_of_genes;
if ($use_func_files > 0){

	######################################
	# for bacteria, use genbank files
	######################################
	open(LS,"ls $prot_dir/*gb |");
	while(my $gb_file = <LS>){
		chomp($gb_file);
		open(GB,$gb_file);
		my $current_cog;
		my $current_gene;
		while(my $l=<GB>){
			if ($l=~/(COG\d+)/){
				#$current_cog = $1;
			}
			if ($l =~/protein_id=\"(.*)\"/ && $current_cog){
				$current_gene = $1;
				$cog_of_genes{$current_gene} .= ",$current_cog";
			}
			if ($l =~/locus_tag=\"(.*)\"/ && $current_cog){
                        	$current_gene = $1;
				$cog_of_genes{$current_gene} .= ",$current_cog";
        	        }
		}
		close(GB);
	}
	close(LS);
}

else{
	######################################
        # for eucaryotes, use GFF files
        ######################################
	open(LS,"ls $prot_dir/*gff |");
	while(my $file = <LS>){
		chomp($file);
		open(GFF,$file);
		my $current_GO;
		my $current_gene;
		while(my $l=<GFF>){
			if ($l =~/ID=([^;]+);.*=(GO:[^;]+);/){
				$current_gene = $1;
				$current_gene =~s/:/_/g;
				if ($l=~/=(GO:\d+[^;]+);/){
					$current_GO = $1;
					$cog_of_genes{$current_gene} .= ",$current_GO";
				}
				if ($l=~/Note=([^;]+);/){
					my $function = $1;
					$function =~s/'//g;
					$functions{$current_gene} = $function;
				}
			}
		}
		close(GFF);
	}
	close(LS);
}

my %letters_of_cog;
my %count_letter;
if (scalar keys(%cog_of_genes) > 1){
	open(WHOG,"$dirname/../COG/whog");
	while(my $l=<WHOG>){
		if ($l=~/\[(\w+)\] (COG\d+) /){
			my $letter = $1;
			my $cogid = $2;
			$letters_of_cog{$cogid} = $letter;
		}
	}
	close(WHOG);

	foreach my $gene(keys(%cog_of_genes)){
		my $cluster = $cluster_of_gene{$gene};
		my $cogs = $cog_of_genes{$gene};
		my $function = $functions{$gene};
		$cogs_of_cluster{$cluster} = $cogs;
		$functions{$cluster} = $function;	
	}
}
if (scalar keys(%cogs_of_cluster) > 1) {
	open(COG_CLUST,">$cog_clusters");	
	foreach my $cluster(sort {$a<=>$b} keys(%cogs_of_cluster)){
		my $cogids = $cogs_of_cluster{$cluster};
		my @ids = split(/,/,$cogids);
		foreach my $id(@ids){
			if ($id =~/\w+/ && $cluster){
				my $letter = $letters_of_cog{$id};
				my @letters = split(//,$letter);
				my $letters_string = join("\t",@letters);
				#$letter =~s//\t/g;
				print COG_CLUST "$cluster	$id	$letters_string\n";
			}
		}
	}
	close(COG_CLUST);
}
else{

	my $is_plus = `which rpsblast+`;
	my $is_not_plus = `which rpsblast`;
	my $software = "rpsblast";
	if ($is_plus){$software = "rpsblast+";}
	system("$software -query $pav_matrix.selection_prot.fa -db $dirname/../COG/Cog -out $pav_matrix.selection.rps-blast.out -evalue 1e-2 -outfmt '7 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs'");

	system("perl $dirname/../COG/bac-genomics-scripts/cdd2cog/cdd2cog.pl -r $pav_matrix.selection.rps-blast.out -c $dirname/../COG/cddid.tbl -f $dirname/../COG/fun.txt -w $dirname/../COG/whog -a");
	system("cp -rf results/protein-id_cog.txt $cog_clusters");
	system("rm -rf results");
}

open(COG,">$cog_outfile");
my %cogs;
open(C,$cog_clusters);
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
close(COG);



my %cogs_of_clusters;
my %cogcats_of_clusters;
open(C,$cog_clusters);
while(<C>){
        my $line = $_;
        $line =~s/\n//g;$line =~s/\r//g;
        my @infos = split(/\t/,$line);
        my $cluster = $infos[0];
        my $cog = $infos[1];
        my $cog_category = $infos[2];
        $cogs_of_clusters{$cluster}{$cog} = 1;
        $cogcats_of_clusters{$cluster} = $cog_category;
}
close(C);

open(CC,">$cog_clusters.2");
foreach my $cluster(sort{$a<=>$b} keys(%genes_of_cluster)){
        my $cog_category = "Unknown";
        my $cog_list = "Unknown";
        if ($cogcats_of_clusters{$cluster}){
                $cog_category = $cogcats_of_clusters{$cluster};
                my $ref_cogs_of_clusters = $cogs_of_clusters{$cluster};
                $cog_list = join(",",keys(%$ref_cogs_of_clusters));
        }
	my $function = "Unknown";
	if ($functions{$cluster}){
		$function = $functions{$cluster};
	}
        if ($accessory_clusters{$cluster}){
                print CC "$cluster\t$cog_list\t$cog_category\t$function\n";
        }
}
close(CC);

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

