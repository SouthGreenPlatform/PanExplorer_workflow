#!/usr/bin/perl

use strict;

my $indir = $ARGV[0];
my $matrix = $ARGV[1];
my $out = $ARGV[2];
my $strain_names = $ARGV[3];

my %strains_of_gb;
open(F,$strain_names);
while(<F>){
	my $line = $_;
	$line =~s/\n//g;$line =~s/\r//g;
	my ($gb,$strain) = split(/\t/,$line);
	$strains_of_gb{$gb} = $strain;
}
close(F);

my @strains;
my %corr;
open(D,"ls $indir/*pep |");
while(<D>){
	my $file = $_;
	my $prot_num = 0;
	my $strain;
	if ($file =~/\/([^\/]*).pep/){
		$strain = $1;
		$strain = $strains_of_gb{$strain};
		push(@strains,$strain);
	}
	open(F,"$file");
        while(<F>){
		if (/>(.*)/){
                        my $prot_id = $1;
                        $prot_num++;
                        my $new_id = "$strain"."_".$prot_num;
			$corr{$new_id} = $prot_id;
                }
        }
        close(F);
}
close(D);

my $cl_num = 0;
my $nb_strains = 1;
open(O,">$out");
open(U,">$out.upsetr.txt");
open(M,">$out.accessory_01matrix.txt");
open(F,$matrix);
print O "ClutserID";
print U "ClutserID";
print M "Gene";
my %hash_place_strains;
my $num_cell = 0;
foreach my $strain(@strains){
        $num_cell++;
        my @words = split(/_/,$strain);
        my $genus = $words[0];
        my $species = $words[1];
        my $shortname = substr($genus,0,3) . "_". substr($species,0,2);
        for (my $j = 2; $j <= $#words; $j++){
                $shortname.="_".$words[$j];
        }
        $shortname = substr($shortname,0,25);
        print O "\t".$strain;
        print U "\t".$shortname;
        print M "\t".$shortname;
        $hash_place_strains{$strain} = $num_cell;
        $nb_strains++;
}
print O "\n";
print U "\n";
print M "\n";
while(<F>){
        $cl_num++;
        my $line = $_;
        $line =~s/\n//g;$line =~s/\r//g;
        my @infos = split(/ /,$line);
        my %cells;
        for (my $i = 1; $i <= $#infos; $i++){
                my $new_id = $infos[$i];
                my $prot_id = $corr{$new_id};
                my $strain;
                if ($new_id =~/^(.*)_\d+$/){$strain=$1;}
                my $num_cell = $hash_place_strains{$strain};
		#print "$strain $num_cell $prot_id $new_id\n";
                $cells{$strain}.= $prot_id.",";
        }
        print O $cl_num;
        print U $cl_num;
        my $concat_accessory = "";
        foreach my $strain(@strains){
                my $val;
                if ($cells{$strain}){
                        $val = $cells{$strain};
                        chop($val);
                }
                else{
                        $val = "-";
                }
                if ($val =~/\w+/){
                        print U "\t1";
                        $concat_accessory .= "\t1";
                }
                else{
                        print U "\t0";
                        $concat_accessory .= "\t0";
                }
                my $concat = $val;
                print O "\t".$concat;
        }
        if ($concat_accessory =~/0/){
                print M $cl_num.$concat_accessory."\n";
        }
        print O "\n";
        print U "\n";
}
close(F);
close(O);
