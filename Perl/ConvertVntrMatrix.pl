#!/usr/bin/perl

use strict;

my $matrix = $ARGV[0];
my $out = $ARGV[1];
my $strain_names = $ARGV[2];
my $bed_dir = $ARGV[3];
my $fasta_dir = $ARGV[4];

my %strains_of_gb;
open(F,$strain_names);
while(<F>){
	my $line = $_;
	$line =~s/\n//g;$line =~s/\r//g;
	my ($gb,$strain) = split(/\t/,$line);
	$strains_of_gb{$gb} = $strain;
}
close(F);

#######################################
# get flanking sequences of repeats
#######################################
my %flankings;
open(BED_DIR," ls $bed_dir/*bed |");
while(<BED_DIR>) {
	my $bedfile = $_;
	if ($bedfile =~/\/([^\/]*)\.fasta/){
		my $accession = $1;
		open(FASTA,"$fasta_dir/$accession.fasta");
		my $id = "";
		my %sequences;
		while(<FASTA>){
			if (/>([^\s]+)/){
				$id = $1;
				if ($id =~/(.*)\.\d+/){$id=$1;}
			}
			else{
				my $line = $_;
				$line =~s/\n//g;$line =~s/\r//g;
				$sequences{$id}.=$line;
			}
		}
		close(FASTA);

		open(BED,$bedfile);
		while(<BED>){
			my $line = $_;
			$line =~s/\n//g;$line =~s/\r//g;
			my ($accession2,$start,$end,$repeat) = split(/\t/,$line);
			my ($repeat_name,$repeat,$nb_repeat) = split(/_/,$repeat);
			my $sequence = $sequences{$accession2};
			my $flanking1 = substr($sequence,$start-300,300);
			my $flanking2 = substr($sequence,$end,300);
			my $strain = $strains_of_gb{$accession};
			my $flanking1reverse = reverse $flanking1;
		        $flanking1reverse =~ tr/ATGCatgc/TACGtacg/;
			my $flanking2reverse = reverse $flanking2;
			$flanking2reverse =~ tr/ATGCatgc/TACGtacg/;
			$flankings{$repeat_name}{$strain}{$flanking1reverse} = 1;
			$flankings{$repeat_name}{$strain}{$flanking2reverse} = 1;
			$flankings{$repeat_name}{$strain}{$flanking1} = 1;
			$flankings{$repeat_name}{$strain}{$flanking2} = 1;
		}
		close(BED);
	}
}
close(BED_DIR);


my $nb_strains = 1;
open(O,">$out");
open(F,$matrix);
my $first_line = <F>;
$first_line =~s/\n//g;$first_line =~s/\r//g;
my @infos_first_line = split(/\t/,$first_line);
$first_line =~s/ClutserID/ID\tRepeat\tFlanking/g;
print O "$first_line\n";
my $n = 0;
my $monomorphic = 0;
my $m = 0;
while(my $line = <F>){
	chomp($line);
	my @infos = split(/\t/,$line);
       	my %hash;
	my %alleles;
	my $concat = "";
	my $concat_flanking = "";
	my %repeats;
	my %flanking_sequences;	
	for (my $i = 1; $i <= $#infos; $i++){
		if ($infos[$i] !~/-/){
			my ($repeat_name,$repeat,$nb_repeat) = split(/_/,$infos[$i]);
			if ($nb_repeat =~/(\d+)\.+/){$nb_repeat = $1;}
			my $size_repeat = length($repeat);
			$concat .= "\t".$nb_repeat;
			$repeats{$repeat}++;
			$alleles{$nb_repeat}++;
			$hash{$size_repeat}++;
			my $genome_name = $infos_first_line[$i];
			
			my $ref_subhash = $flankings{$repeat_name}{$genome_name};
			my %subhash = %$ref_subhash;
			foreach my $flank(keys(%subhash)){
				$flanking_sequences{$flank} = 1;
			}
		}
		else{
			$concat .= "\t".$infos[$i];
		}
	}
	my $concat_flanking = join(",",keys(%flanking_sequences));
	foreach my $flank(keys(%flanking_sequences)){
		
	}
	# different patterns
	if (scalar keys(%hash) > 1){
		#print O join(",",keys(%repeats)).":".$concat_flanking."$concat\n";
		$m++;
	}
	# monomorphic
	elsif (scalar keys(%alleles) == 1){
		#print O $infos[0]. "\t". join(",",keys(%repeats))."\t".$concat_flanking."$concat\n";
		$monomorphic++;
	}
	# found several times for the same sample
	elsif ($line =~/,/){
		#print O join(",",keys(%repeats)).":".$concat_flanking."$concat\n";	
	}
	else{
		print O $infos[0]. "\t". join(",",keys(%repeats))."\t".$concat_flanking."$concat\n";
		$n++;
	}
}
close(F);
close(O);

