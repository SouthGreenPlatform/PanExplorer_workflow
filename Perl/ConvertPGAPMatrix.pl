#!/usr/bin/perl

use strict;
use FindBin;
use lib $FindBin::Bin;
use GeneId qw(parse_gff_attributes);

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
                if ($infos[2] eq 'CDS'){
                        # map the tool identifier (GFF ID) to the gene identifier used everywhere
                        # else in PanExplorer: locus_tag if available, protein_id otherwise (see GeneId.pm)
                        my $attr = parse_gff_attributes($infos[8]);
                        my $id = $attr->{'ID'};
                        my $gene = defined $attr->{'locus_tag'} ? $attr->{'locus_tag'} : $attr->{'protein_id'};
                        if (defined $id && defined $gene){
                                $corr{$id} = $gene;
                        }
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
