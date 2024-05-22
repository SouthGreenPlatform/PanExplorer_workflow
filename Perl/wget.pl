#!/usr/bin/perl

use strict;

use File::Basename;
my $dirname = dirname(__FILE__);

system("wget https://ftp.ncbi.nlm.nih.gov/genomes/GENOME_REPORTS/prokaryotes.txt");
system("wget https://ftp.ncbi.nlm.nih.gov/genomes/GENOME_REPORTS/eukaryotes.txt");

my %continents;
open(F,"countries.txt");
<F>;
while(my $line =<F>){
        chomp($line);
        my ($continent,$country) = split(/,/,$line);
        $continents{$country} = $continent;
}
close(F);

my $input = $ARGV[0];
my $outdir = $ARGV[1];
my $private_genomes = $ARGV[2];

system("cat $input $private_genomes >$outdir/list.txt");

my $concat = "";
open(O2,">$outdir/genbanks.txt");
open(O,">$outdir/strains.txt");
open(GENES,">$outdir/genes.txt");
open(L,">$outdir/list_genomes.txt");
open(L2,">$outdir/list_genomes2.txt");
open(L3,">$outdir/genomes.txt");
open(L4,">$outdir/genomes2.txt");
open(SEQFILE,">$outdir/seqfile");
open(PanSN,">$outdir/all_genomes.fa");
open(METADATA,">$outdir/metadata_strains.txt");

open(F,"$outdir/list.txt");
#open(TEST,">$outdir/test");
while(my $line =<F>){
	chomp($line);
	my $genbank = $line;
	if (!-e "$genbank"){
		my $grep = `grep '$line' prokaryotes.txt`;
		#print "$genbank $line aaa $grep\n";exit;
		my @infos = split(/\t/,$grep);
	        my $status = $infos[15];
        	if ($status !~/Complete Genome/ && $status !~/Chromosome/){
	                #next;
        	}
	        my $ftp_path = $infos[$#infos -2];

        	$ftp_path =~s/ftp:/http:/g;
	        my @table = split(/\//,$ftp_path);
        	my $name = $table[$#table];
	        my $prot_file = "$ftp_path/$name"."_protein.faa.gz";
        	my $gbff = "$ftp_path/$name"."_genomic.gbff.gz";
	        my $gff = "$ftp_path/$name"."_genomic.gff.gz";
        	my $genome_fasta = "$ftp_path/$name"."_genomic.fna.gz";
	        my @particules = split(/_/,$name);

        	`wget -O $outdir/$genbank.fasta.gz $genome_fasta`;
		`gunzip $outdir/$genbank.fasta.gz`;
		`wget -O $outdir/$genbank.gb.gz $gbff`;
	        system("gunzip $outdir/$genbank.gb.gz");


		################################################################
		# for eukaryotes
		################################################################
		if (!$grep){
			$grep = `grep '$line' eukaryotes.txt`;
			my @infos = split(/\t/,$grep);
			my $gca = $infos[8];
			if ($gca =~/GCA_(\d\d\d)(\d\d\d)(\d\d\d)/){
				my $part1 = $1;
				$ftp_path = "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/$1/$2/$3";
				`wget -O $outdir/$gca.index.html $ftp_path`;
				my $grep_name = `grep '$gca' $outdir/$gca.index.html`;
				unlink("$outdir/$gca.index.html");
				my $name;
				if ($grep_name =~/href=\"(.*)\"/){
					$name = $1;
					$name=~s/\///g;
				}
				$ftp_path = $ftp_path."/$name";
				my $prot_file = "$ftp_path/$name"."_protein.faa.gz";
				my $gbff = "$ftp_path/$name"."_genomic.gbff.gz";
				my $gff = "$ftp_path/$name"."_genomic.gff.gz";
				my $genome_fasta = "$ftp_path/$name"."_genomic.fna.gz";
				my @particules = split(/_/,$name);

                `wget -O $outdir/$genbank.fasta.gz $genome_fasta`;
                `gunzip $outdir/$genbank.fasta.gz`;
                `wget -O $outdir/$genbank.gb.gz $gbff`;
		`wget -O $outdir/$genbank.faa.gz $prot_file`;
                system("gunzip $outdir/$genbank.gb.gz");
		system("gunzip $outdir/$genbank.faa.gz");

			}
		}
	}
	else{
		my $genbank_file = $genbank;
		my $grep = `grep 'LOCUS' $genbank_file`;
		$genbank = "unknown";
		if ($grep =~/LOCUS\s+([\-\:\w]+)/){$genbank = $1;}

		#$genbank =~s/\:/_/g;	

		my $cmd = "cp -rf $genbank_file $outdir/$genbank.gb";
		system($cmd);

		my %genome_seqs;
		my $current_chr;
		my $go = 0;
		open(G,"$outdir/$genbank.gb");
		while(<G>){
			if ($go == 1 && /(\d+) (.*)$/){
				my $line = $2;
				$line =~s/ //g;
				$genome_seqs{$current_chr}.=$line;
			}
			if (/LOCUS       ([^\s]+)/){
				$current_chr = $1;
			}
			if (/ORIGIN/){$go = 1;}
			if (/^\/\//){$go = 0;}
		}
		close(G);

		open(FASTA,">$outdir/$genbank.fasta");
		foreach my $ch(keys(%genome_seqs)){
			print FASTA ">$ch\n";
			my $seq = $genome_seqs{$ch};
			print FASTA "$seq\n";
		}
		close(FASTA);
	}
	#my $get_organism_line = `head -10 $outdir/$genbank.gb | grep DEFINITION `;
	my $get_organism_line = `head -10 $outdir/$genbank.gb | grep -A 1 DEFINITION `;

	# if several lines for DEFINITION, concatenate the lines
	my @lines_organism = split(/\n/,$get_organism_line);
	my $first_line = $lines_organism[0];
	my $second_line = $lines_organism[1];
	if ($second_line =~/^            (.*)/){
		$get_organism_line = $first_line. " ".$1;
	}
	else{
		$get_organism_line = $first_line;
	}

        my $strain;
	my $genus;
        if ($get_organism_line =~/DEFINITION  (.*)$/){
                $strain = $1;
		($genus) = split(/\s/,$strain);
        }
	my $country = `head -100 $outdir/$genbank.gb | grep country `;
        $country =~s/^\s+//g;
        $country =~s/\/country=//g;
        $country =~s/\"//g;
        $country =~s/\n//g;$country =~s/\r//g;
        if ($country =~/:/){
                my $city;
                ($country,$city) = split(/:/,$country);
        }
        if ($country eq ""){$country = "unresolved";}
        my $continent = "unresolved";
        if ($continents{$country}){
                $continent = $continents{$country};
        }
	$strain =~s/\.//g;
	my ($info1,$info2 ) = split(",",$strain);
        $strain = $info1;
        $strain =~s/ /_/g;
        $strain =~s/strain_//g;
        $strain =~s/_chromosome//g;
        $strain =~s/_genome//g;
        $strain =~s/_str_/_/g;
	$strain =~s/[^\w\-\_]//g;
	$strain =~s/\-/_/g;

	print O "$genbank	$strain\n";	
	$concat .= "$genbank,";
	print L "$genbank	$outdir/$genbank.gb\n";
	print L2 "$genbank\n";
	print L3 "$outdir/$genbank.fasta\n";
	print L4 "$outdir/$genbank.fasta	$strain\n";
	print SEQFILE "$genbank	$outdir/$genbank.fasta\n";
	print METADATA "$strain\t$genus\t$country\t$continent\n";


	my %genome_sequences;
	my $genome = "";
	my $seqid = "";
        open(GENOME,"$outdir/$genbank.fasta");
        while(<GENOME>){
                if (!/^>/){
                        my $line = $_;
                        $line =~s/\n//g;$line =~s/\r//g;
			$genome_sequences{$seqid} .= $line;
                        $genome .= $line;
			print PanSN $_;
                }
		elsif (/>([^\s]+)/){
			$seqid = $1;
			$seqid =~s/\.\d+$//g;
			print PanSN ">$strain#$seqid\n";
		}
        }
        close(GENOME);


	open(N,">$outdir/$genbank.nuc");
        open(P,">$outdir/$genbank.pep");
        open(FUNC,">$outdir/$genbank.func");
        my $go = 0;
        my $start;
        my $end;
        my $product;
        my $complement = 0;
        my $end_gene = "no";
        my $protein = "";

	#`sed -i "s/'//g" $outdir/$genbank.gb`;

	my $has_translation = `grep -c 'translation=' $outdir/$genbank.gb`;
	$has_translation =~s/\n//g;$has_translation =~s/\r//g;
        open(G,"$outdir/$genbank.gb");
        my $current_gene;
	my $current_chrom;
        while(<G>){
                if (/^\s+ORGANISM\s+(.*)$/){
                }
		if (/protein_id=\"(.*)\"/){
			$current_gene = $1;
		}
		if (/LOCUS       ([^\s]+)/){
			$current_chrom = $1;
		}        
                if (/locus_tag=\"(.*)\"/){
                        $current_gene = $1;
		}
                if ($go == 1){
                        my $line = $_;
                        $line =~s/ //g;
                        $line =~s/\n//g;$line =~s/\r//g;
                        $protein .= $line;
                        if (/\"$/){
                                $protein =~s/\"//g;
                                $end_gene = "yes";

                        }
                }
                if (/\/translation=\"(.*)/ or ($has_translation == 0 && /protein_id=/)){
                        $go = 1;
                        $protein .= $1;
			print P ">$current_gene\n";
			print N ">$current_gene\n";
			print GENES "$current_gene $product [$strain]\n";

                        if ($protein =~/\"$/){
                                $end_gene = "yes";
                        }
			if ($has_translation == 0){$end_gene = "yes";}
                }
                if ($end_gene eq "yes"){
                        $protein =~s/\"//g;
                        print P "$protein\n";
                        $protein = "";
                                my $length = $end - $start + 1;

				#print "okkkk $current_chrom\n";
                                #my $geneseq = substr($genome,$start-1,$length);
                                my $geneseq = substr($genome_sequences{$current_chrom},$start-1,$length);


                                if ($complement){
                                        my $revcomp = reverse $geneseq;
                                        $revcomp =~ tr/ATGCatgc/TACGtacg/;
                                        $geneseq = $revcomp;
                                }
                                print N "$geneseq\n";
                                print FUNC "$current_gene  -       $product\n";
                                $go = 0;
                        $end_gene = "no";
                }
                if (/\/product=\"(.*)\"/){
                        $product = $1;
                }
                if (/^\s+CDS\s+(\d+)\.\.(\d+)\s*$/){
                        $start = $1;
                        $end = $2;
                        $complement = 0;
                }
                if (/^\s+CDS\s+complement\((\d+)\.\.(\d+)\)\s*$/){
                        $start = $1;
                        $end = $2;
                        $complement = 1;
                }
        }
        close(G);
        close(P);
        close(N);
        close(FUNC);

	if ($has_translation == 0){
		system("perl $dirname/translate.pl $outdir/$genbank.nuc $outdir/$genbank.pep");
	}
	
	my $prot_num = 0;
        open(PRT,">$outdir/$genbank.prt");
        open(P,"$outdir/$genbank.pep");
        while(<P>){
                if (/>(.*)/){
                        my $prot_id = $1;
                        $prot_num++;
                        my $new_id = "$strain"."_".$prot_num;
                        print PRT ">$new_id\n";
                }
                else{
                        print PRT $_;
                }
        }
        close(P);
        close(PRT);
}
close(F);
close(O);
close(METADATA);
chop($concat);
print O2 $concat;
close(O2);
close(L);
close(L2);
close(L3);
close(L4);
close(GENES);
close(SEQFILE);
close(PanSN);
#close(TEST);

unlink("prokaryotes.txt");
unlink("eukaryotes.txt");
