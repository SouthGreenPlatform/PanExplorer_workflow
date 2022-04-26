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
open(F,$matrix);
my $firstline = <F>;
$firstline =~s/\n//g;$firstline =~s/\r//g;
my @infos = split(/\t/,$firstline);
print O "ClutserID";
print U "ClutserID";
print M "Gene";
for (my $j=1; $j <= $#infos; $j++){
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
while(<F>){
	print O $_;
}
close(F);
close(O);
