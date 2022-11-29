#!/usr/bin/perl

use strict;

my $nb_file_found = 0;
open(LS,"ls *Nfowleri.snp |");
while(<LS>){
	$nb_file_found++;
}
close(LS);
if ($nb_file_found < 1){
	print "No file found with extension \".Nfowleri.snp\". This script allows to merge multiple SNP files (named as follows \"strainname.Nfowleri.snp\") obtained with VarScan into a global VCF file.\n";
	exit;
}


print "##fileformat=VCFv4.1\n";
print "#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	";

my %corr;
open(F,"/media/results_datacalcul/Alexis/Naegleria/raw/SRR.info.xls");
while(<F>){
        my $line = $_;
        $line =~s/\n//g;$line =~s/\r//g;
        my ($srr,$name) = split(/\t/,$line);
        #$corr{$srr} = $name;
}
close(F);

my %str;
my @strains;
my %hash;
my %genotypes;
my %alternates;
my %reference_bases;
open(LS,"ls *Nfowleri.snp |");
#open(LS,"ls *Nlovaniensis.snp |");
#open(LS,"ls ../../tuberculosis/VCF_comparaison_avec_reference/*vcf |");
while(<LS>){
	my $file = $_;
	$file =~s/\n//g;$file =~s/\r//g;
	my $strain;
	if ($file =~/^(.*).Nfowleri.snp/){
	#if ($file =~/^(.*).Nlovaniensis.snp/){
	#if ($file =~/\/([^\/]*)\.vcf/){
		$strain = $1;
		$strain =~s/_/-/g;
		#if ($strain =~/SRR12781212/ or $strain =~/SRR12781213/ or $strain =~/SRR12781214/ or $strain =~/SRR12781215/ or $strain =~/SRR12781216/){next;}
		#if ($strain !~/NF/){next;}
		push(@strains,$strain);
		$str{$strain} = 1;
		#open(F,"$strain.Klebsiella_pneumoniae_subsp._pneumoniae_HS11286.snp");
		open(F,$file);
		<F>;
		while(<F>){
			my $line = $_;
			$line =~s/\n//g;$line =~s/\r//g;
			if (/^#/){next;}
			my ($chr,$pos,$ref,$alt,$info) = split("\t",$line);
			my ($Cons,$Cov,$Reads1,$Reads2,$Freq) = split(":",$info);
 			$Freq =~s/%//g;$Freq =~s/,/./g;
			$chr =~s/\.1//g;
			$alternates{$chr}{$pos}{$alt}=1;
			$reference_bases{$chr}{$pos} = $ref;
			if ($Cov > 6 && $Freq > 70){
				$hash{$chr}{$pos}{$strain} = "1/1";
			}
			elsif ($Cov > 6 && $Freq > 34){
				$hash{$chr}{$pos}{$strain} = "0/1";
			}
			else{
                                #$hash{$chr}{$pos}{$strain} = "0/0";
                        }
		}
		close(F);
	}
}
close(LS);

my %depths;
open(LS,"ls *.depth |");
while(<LS>){
	my $file = $_;
	$file =~s/\n//g;$file =~s/\r//g;
	my $strain;
	if ($file =~/^(.*).Nfowleri.depth/){
	#if ($file =~/^(.*).Nlovaniensis.depth/){
		$strain = $1;
		$strain =~s/_/-/g;
		if (!$str{$strain}){next;}
		#if ($strain !~/NF/){next;}
		open(F,$file);
		while(<F>){
			my $line = $_;
			$line =~s/\n//g;$line =~s/\r//g;
			my ($chr,$pos,$depth) = split("\t",$line);
			$chr =~s/\.1//g;
			if ($hash{$chr}{$pos}){
				$depths{$chr}{$pos}{$strain}=1;
			}
		}
		close(F);
	}
}
close(LS);

foreach my $strain(@strains){
	if ($corr{$strain}){
		$strain= $corr{$strain};
	}
	print "$strain\t";
}
print "Reference\n";
foreach my $chr(sort {$a<=>$b} keys(%hash)){
my $refsubhash = $hash{$chr};
my %subhash2 = %$refsubhash;
foreach my $pos(sort {$a<=>$b} keys(%subhash2)){
	my $refbase = $reference_bases{$chr}{$pos};
	my $ref_alt = $alternates{$chr}{$pos};
	my %subhash_alt = %$ref_alt;
	my $alternate;

	# discard multiallelic
	if (scalar keys(%subhash_alt) > 1){
		next;
	}
	else{
		foreach my $alt(keys(%subhash_alt)){
			$alternate = $alt;
		}
	}
	
	print "$chr	$pos	.	$refbase	$alternate	186.96	.	PASS	GT:AD:DP:GQ:PL";
	foreach my $strain(@strains){
		my $val = $hash{$chr}{$pos}{$strain};
		if ($val){
			print "	$val";
		}
		else{
			if ($depths{$chr}{$pos}{$strain}){
				print "	0/0";
			}
			else{
				print "	./.";
			}
		}
		
	}
	print "	0/0\n";
}
}


