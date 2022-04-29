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


my %corr;
open(D,"ls $indir/*rmdup.gff |");
while(<D>){
	my $file = $_;
	open(F,"$file");
        while(<F>){
                my @infos = split(/\t/,$_);
                if ($infos[2] eq 'CDS' && /ID=([^;]*);.*protein_id=([^;]*);/){
                        my $id = $1;
                        my $protid = $2;
                        $corr{$id} = $protid;
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
my $firstline = <F>;
$firstline =~s/\n//g;$firstline =~s/\r//g;
my @infos = split(/,/,$firstline);
print O "ClutserID";
print U "ClutserID";
print M "Gene";
for (my $j=14; $j <= $#infos; $j++){
        my $gbfile = $infos[$j];
        $gbfile =~s/\"//g;
        $gbfile =~s/\.gb\.filt//g;
	$gbfile =~s/\.gb\.rmdup//g;
	
        my $strain = $strains_of_gb{$gbfile};
        print O "\t".$strain;
        print U "\t".$strain;
        print M "\t".$strain;
        $nb_strains++;
}
print O "\n";
print U "\n";
print M "\n";
while(<F>){
        $cl_num++;
        my $line = $_;
        $line =~s/\n//g;$line =~s/\r//g;
        my @infos = split(/,/,$line);
        print O $cl_num;
        print U $cl_num;
        my $concat_accessory = "";
        for (my $i = 14; $i <= $#infos; $i++){
                my $val = $infos[$i];
                $val =~s/\"//g;
                if ($val =~/\w+/){
                        print U "\t1";
                        $concat_accessory .= "\t1";
                }
                else{
                        print U "\t0";
                        $concat_accessory .= "\t0";
                }
                my @genes = split(/;/,$val);
                my $concat = "";
                foreach my $gene(@genes){
                        my $prot_id = $corr{$gene};
			if (!$prot_id){$prot_id = $gene;}
                        $concat .= "$prot_id,"
                }
                chop($concat);
                if (scalar @genes == 0){
                        $concat = "-";
                }
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
close(U);
close(M);
