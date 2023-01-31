#!/usr/bin/perl

use strict;

system("wget https://ftp.ncbi.nlm.nih.gov/genomes/GENOME_REPORTS/prokaryotes.txt");

my %continents;
open(F,"countries.txt");
<F>;
while(<F>){
        my $line = $_;
        $line =~s/\n//g;$line =~s/\r//g;
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
open(METADATA,">$outdir/metadata_strains.txt");
open(F,"$outdir/list.txt");
while(<F>){
	my $line = $_;$line =~s/\n//g;$line =~s/\r//g;
	my $genbank = $line;
	
	if (!-e "$genbank"){
		my $grep = `grep $line prokaryotes.txt`;
		my @infos = split(/\t/,$grep);
	        my $status = $infos[15];
        	if ($status !~/Complete Genome/ && $status !~/Chromosome/){
	                next;
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
	}
	else{
		my $genbank_file = $genbank;
		my $grep = `grep 'LOCUS' $genbank_file`;
		$genbank = "unknown";
		if ($grep =~/LOCUS\s+(\w+)/){$genbank = $1;}

		my $cmd = "cp -rf $genbank_file $outdir/$genbank.gb";
		system($cmd);

		my $go = 0;
		open(FASTA,">$outdir/$genbank.fasta");
		print FASTA ">$genbank\n";
		open(G,"$outdir/$genbank.gb");
		while(<G>){
			if ($go == 1 && /(\d+) (.*)$/){
				my $line = $2;
				$line =~s/ //g;
				print FASTA $line;
			}
			if (/ORIGIN/){$go = 1;}
		}
		close(G);
		close(FASTA);
	}
	
	my $get_organism_line = `head -10 $outdir/$genbank.gb | grep DEFINITION `;
        my $strain;
	my $genus;
        if ($get_organism_line =~/DEFINITION  (.*)$/){
                $strain = $1;
		($genus) = split(/\s/,$strain);
        }
	my $country = `grep country $outdir/$genbank.gb`;
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
	my ($info1,$info2 ) = split(",",$strain);
        $strain = $info1;
        $strain =~s/ /_/g;
        $strain =~s/strain_//g;
        $strain =~s/_chromosome//g;
        $strain =~s/_genome//g;
        $strain =~s/str\._//g;
	$strain =~s/\=//g;
	$strain =~s/\///g;
        $strain =~s/\(//g;
        $strain =~s/\)//g;
        $strain =~s/\.//g;
        $strain =~s/\-/_/g;
	$strain =~s/\:/_/g;
	print O "$genbank	$strain\n";	
	$concat .= "$genbank,";
	print L "$genbank	$outdir/$genbank.gb\n";
	print L2 "$genbank\n";
	
	print METADATA "$strain\t$genus\t$country\t$continent\n";


	my $genome = "";
        open(GENOME,"$outdir/$genbank.fasta");
        while(<GENOME>){
                if (!/^>/){
                        my $line = $_;
                        $line =~s/\n//g;$line =~s/\r//g;
                        $genome .= $line;
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
        open(G,"$outdir/$genbank.gb");
        my $current_gene;
        while(<G>){
                if (/^\s+ORGANISM\s+(.*)$/){
                }
		if (/protein_id=\"(.*)\"/){
			$current_gene = $1;
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
                if (/\/translation=\"(.*)/){
                        $go = 1;
                        $protein .= $1;
			print P ">$current_gene\n";
			print N ">$current_gene\n";
			print GENES "$current_gene $product [$strain]\n";

                        if ($protein =~/\"$/){
                                $end_gene = "yes";
                        }

                }
                if ($end_gene eq "yes"){
                        $protein =~s/\"//g;
                        print P "$protein\n";
                        $protein = "";
                                my $length = $end - $start + 1;
                                my $geneseq = substr($genome,$start-1,$length);


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
                if (/^\s+CDS\s+(\d+)\.\.(\d+)$/){
                        $start = $1;
                        $end = $2;
                        $complement = 0;
                }
                if (/^\s+CDS\s+complement\((\d+)\.\.(\d+)\)$/){
                        $start = $1;
                        $end = $2;
                        $complement = 1;
                }
        }
        close(G);
        close(P);
        close(N);
        close(FUNC);
	
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
close(GENES);

unlink("prokaryotes.txt");
