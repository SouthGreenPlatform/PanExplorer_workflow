#!/usr/bin/perl

use strict;
use File::Basename;
my $dirname = dirname(__FILE__);

my $pav_matrix = $ARGV[0];
my $prot_dir = $ARGV[1];
my $cog_outfile = $ARGV[2];
my $kegg = $ARGV[3];
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
			$id =~s/://g;
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
my %keggs_of_cluster;
my %functions;
my %accessory_clusters;
my %genes_of_cluster;
my %cluster_of_gene;
open(S,">$pav_matrix.selection_prot.fa");
open(O,$pav_matrix);
<O>;
my $n = 0;
my %selection_protein;
while(<O>){
	$n++;
        my $line = $_;
        $line =~s/\n//g;$line =~s/\r//g;
        my @infos = split(/\t/,$line);
        my $cluster = $infos[0];
	#print S ">$cluster\n";
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
				#print S $proteins{$genename};
				$selection_protein{$cluster} = $proteins{$genename};
				print S ">$cluster\n";print S $proteins{$genename}."\n";
				$cluster_of_gene{$genename} = $cluster;
				#if ($use_func_files > 0){
				#	my $function = `grep '$genename' $prot_dir/*func`;
				#        $function =~s/\n//g;$function =~s/\r//g;
				#        my @t = split(/$genename/,$function);
				#        if ($t[1] =~/-\s+(.*)/){
				# 	        $function = $1;
				#        	$function =~s/'//g;
				#	}
				#        $functions{$cluster} = $function;
				#	last LINE;
				#}
			}
		}
        }
}
close(O);


#########################################################################################################################################
# Test if COG are already provided in Genbank files. If provided, rpsblast is not launched, we get the COG information from Genbank files
#########################################################################################################################################
my %cog_of_genes;
my %kegg_of_genes;
if ($use_func_files > 0){
	######################################
	# for bacteria, use genbank files
	######################################
	open(LS,"ls $prot_dir/*gb |");
	while(my $gb_file = <LS>){
		chomp($gb_file);
		open(GB,$gb_file);
		my $current_cog;
		my $current_kegg;
		my $current_gene;
		while(my $l=<GB>){
			if ($l=~/(COG\d+)/){
				$current_cog = $1;
				if ($current_gene){
					$cog_of_genes{$current_gene}{$current_cog} = 1;
				}
			}
			if ($l =~/     CDS             /){
				$current_cog = "";
				$current_gene = "";
				$current_kegg = "";
			}
			if ($l =~/protein_id=\"(.*)\"/){
				$current_gene = $1;
				$current_gene =~s/://g;
				if ($current_cog){
					$cog_of_genes{$current_gene}{$current_cog} = 1;
				}
				if ($current_kegg){
                                        $kegg_of_genes{$current_gene}{$current_kegg} = 1;
                                }
			}
			if ($l =~/KEGG:(K\d+)/){
				$current_kegg = $1;
				if ($current_gene){
					$kegg_of_genes{$current_gene}{$current_kegg} = 1;
				}
			}
			if ($l =~/locus_tag=\"(.*)\"/){
                        	$current_gene = $1;
				$current_gene =~s/://g;
				if ($current_cog){
					$cog_of_genes{$current_gene}{$current_cog} = 1;
				}
				if ($current_kegg){
					$kegg_of_genes{$current_gene}{$current_kegg} = 1;
				}
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
			$l =~s/\n//g;$l =~s/\r//g;
			if ($l =~/ID=([^;]+);.*=(GO:[^;]+);*/){
				$current_gene = $1;
				$current_gene =~s/:/_/g;
				if ($l=~/=(GO:\d+[^;]+);*/){
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
		my @cogs;
		if ($cog_of_genes{$gene}){
			my $ref_hash = $cog_of_genes{$gene};
			my %subhash = %$ref_hash;
			@cogs = keys(%subhash);
		}
		
		my $function = $functions{$gene};
		foreach my $cog(@cogs){
			$cogs_of_cluster{$cluster}{$cog}++;
		}
		$functions{$cluster} = $function;	
	}
}
if (scalar keys(%kegg_of_genes) > 1){
	foreach my $gene(keys(%kegg_of_genes)){
		my $cluster = $cluster_of_gene{$gene};
                my @keggs;
                if ($kegg_of_genes{$gene}){
                        my $ref_hash = $kegg_of_genes{$gene};
                        my %subhash = %$ref_hash;
                        @keggs = keys(%subhash);
                }

                foreach my $kegg(@keggs){
                        $keggs_of_cluster{$cluster}{$kegg}++;
                }
	}
}

open(KEGG,">$kegg");
if (scalar keys(%keggs_of_cluster) > 1) {

	my $DIR_KEGG = "/usr/local/bin/PanExplorer_workflow/KEGG";
	my %pathways;
	my %pathways_of_ko;
	open(KO_TO_PATHWAY,"$DIR_KEGG/ko_pathways.txt");
	while(<KO_TO_PATHWAY>){
		my $line = $_;
		$line =~s/\n//g;$line =~s/\r//g;
		if ($line =~/path:ko/){next;}
		my ($ko,$pathway) = split("\t",$line);
		my ($ko,$koid) = split(":",$ko);
		my ($map,$mapid) = split(":",$pathway);
		$pathways{$mapid}{$koid} = 1;
		$pathways_of_ko{$koid}{$mapid} = 1;
	}
	close(KO_TO_PATHWAY);

	my %pathway_names;
	open(PATHWAYS,"$DIR_KEGG/pathways.txt");
	while(<PATHWAYS>){
		my $line = $_;
		$line =~s/\n//g;$line =~s/\r//g;
		my ($mapid,$pathway_name) = split("\t",$line);
		$pathway_names{$mapid} = $pathway_name;
	}
	close(PATHWAYS);

	open(KO_TO_MODULE,"$DIR_KEGG/ko_modules.txt");
        while(<KO_TO_MODULE>){
                my $line = $_;
                $line =~s/\n//g;$line =~s/\r//g;
                my ($ko,$module) = split("\t",$line);
                my ($ko,$koid) = split(":",$ko);
                my ($module,$moduleid) = split(":",$module);
                $pathways{$moduleid}{$koid} = 1;
                $pathways_of_ko{$koid}{$moduleid} = 1;
        }
        close(KO_TO_PATHWAY);

        my %module_names;
        open(MODULES,"$DIR_KEGG/modules.txt");
        while(<MODULES>){
                my $line = $_;
                $line =~s/\n//g;$line =~s/\r//g;
                my ($id,$module_name) = split("\t",$line);
                $pathway_names{$id} = $module_name;
        }
        close(MODULES);

	my %kegg_infos;
	foreach my $cluster(sort {$a<=>$b} keys(%genes_of_cluster)){
		if ($keggs_of_cluster{$cluster}){
			my $ref_hash = $keggs_of_cluster{$cluster};
			my %subhash = %$ref_hash;
			my @ids = keys(%subhash);
			foreach my $koid(@ids){
				
				if ($pathways_of_ko{$koid}){
					my $ref_hash = $pathways_of_ko{$koid};
					my %subhash = %$ref_hash;
					foreach my $mapid(keys(%subhash)){
						my $pathway = $pathway_names{$mapid};
						$kegg_infos{$mapid}{$koid} = $cluster;
					}
				}
				else{
					my $mapid = "No mapid";
					$pathway_names{$mapid} = "No associated pathway";
					$pathways{$mapid}{$koid} = 1;
					$kegg_infos{$mapid}{$koid} = $cluster;
				}	
			}
		}
	}

	foreach my $mapid(keys(%kegg_infos)){

		#if ($mapid =~/^ko/){next;}
		my $ref_hash = $pathways{$mapid};
		my %subhash = %$ref_hash;
		foreach my $koid(keys(%subhash)){
			my $pathway_name = $pathway_names{$mapid};
			my $cluster = $kegg_infos{$mapid}{$koid};
			print KEGG "$mapid	$pathway_name	$koid	$cluster\n";
		}
	}

}
close(KEGG);

if (scalar keys(%cogs_of_cluster) > 1) {
	open(COG_CLUST,">$cog_clusters");	
	foreach my $cluster(sort {$a<=>$b} keys(%genes_of_cluster)){
		if ($cogs_of_cluster{$cluster}){
		#my $cogids = $cogs_of_cluster{$cluster};
		my $ref_hash = $cogs_of_cluster{$cluster};
		my %subhash = %$ref_hash;
		my @ids = keys(%subhash);

		my @genes = split(",",$genes_of_cluster{$cluster});
		#print "$cluster ".scalar @ids."\n";
		#my @ids = split(/,/,$cogids);
		foreach my $id(@ids){

			if ($id =~/\w+/ && $cluster =~/\w+/){
				my $letter = $letters_of_cog{$id};
				my @letters = split(//,$letter);
				my $letters_string = join("\t",@letters);
				my $nb_times = $cogs_of_cluster{$cluster}{$id} . "/" . scalar @genes;
				#$letter =~s//\t/g;
				print COG_CLUST "$cluster	$id	$letters_string\n";
			}
		}
		}
		else{
			my $protein = $selection_protein{$cluster};
			print S ">$cluster\n";
			print S "$protein\n";
		}
	}
	close(COG_CLUST);
}
close(S);


# if COG not found
if (scalar keys(%cogs_of_cluster) < 1) {
	my $is_plus = `which rpsblast+`;
	my $is_not_plus = `which rpsblast`;
	my $software = "rpsblast";
	if ($is_plus){$software = "rpsblast+";}
	system("$software -query $pav_matrix.selection_prot.fa -db $dirname/../COG/Cog -out $pav_matrix.selection.rps-blast.out -evalue 1e-2 -outfmt '7 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs'");
	
	
	system("perl $dirname/../COG/bac-genomics-scripts/cdd2cog/cdd2cog.pl -r $pav_matrix.selection.rps-blast.out -c $dirname/../COG/cddid.tbl -f $dirname/../COG/fun.txt -w $dirname/../COG/whog -a");
	system("cat results/protein-id_cog.txt >>$cog_clusters");
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
my %c                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          