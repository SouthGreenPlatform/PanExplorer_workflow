#!/usr/bin/perl

my $in = $ARGV[0];
my $out = $ARGV[1];
my $rowmetadata = $ARGV[2];
my $colmetadata = $ARGV[3];

my %hash;
my %hash2;
my $n = 0;
open(M,$rowmetadata);
while(<M>){
	my $line = $_;
	$line =~s/\n//g;$line =~s/\r//g;
	$n++;
	my @infos = split(/\t/,$line);
	my $genus = $infos[1];
	my $country = $infos[2];
	my $continent = $infos[3];
	$hash{$n} = "Genus: $genus <br>Country: $country <br> Continent: $continent";
}
close(M);

my $k = 0;
open(M,$colmetadata);
while(<M>){
        my $line = $_;
        $line =~s/\n//g;$line =~s/\r//g;
        $k++;
        my @infos = split(/\t/,$line);
        my $cog = $infos[1];
	my $cog_cat = $infos[2];
	my $function = $infos[3];
        $hash2{$k} = "Description: $function <br>COG: $cog <br>COG category: $cog_cat";
}
close(M);

open(T,">$in.metadata_tooltips.txt");
foreach my $n(sort {$a<=>$b} keys(%hash)){
	
	my $country = $hash{$n};
	foreach my $k(sort {$a<=>$b} keys(%hash2)){
		my $cog = $hash2{$k};
		print T "<br>$cog <br><br>$country \t";
	}
	print T "\n";
}
close(T);

use strict;
system("ulimit -s 163840;Rscript $dirname/../R/heatmaply.R -f $in -o $out -c $colmetadata -r $rowmetadata");
