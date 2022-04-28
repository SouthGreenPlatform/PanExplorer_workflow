#!/usr/bin/perl

use strict;

my $infile = $ARGV[0];
my $heatmap = $ARGV[1];

my %samples;
my %hash;
my $cl_num = 0;
my $nb_strains = 1;
open(M,">$heatmap.accessory_01matrix.txt");
open(F,"$infile");
my $firstline = <F>;
$firstline =~s/\n//g;$firstline =~s/\r//g;
#my @infos = split(/","/,$firstline);
my @infos = split(/\t/,$firstline);
print M "Gene";
#for (my $j=14; $j <= $#infos; $j++){
for (my $j=1; $j <= $#infos; $j++){
        my $gbfile = $infos[$j];
        $gbfile =~s/\"//g;
        $gbfile =~s/\.gb\.filt//g;
        print M "\t".$gbfile;
	$samples{$j} = $gbfile;
        $nb_strains++;
}
print M "\n";
while(<F>){
        $cl_num++;
        my $line = $_;
        $line =~s/\n//g;$line =~s/\r//g;
        #my @infos = split(/","/,$line);
        my @infos = split(/\t/,$line);
        my $concat_accessory = "";
        #for (my $i = 14; $i <= $#infos; $i++){
	for (my $i = 1; $i <= $#infos; $i++){
                my $val = $infos[$i];
		my $sample = $samples{$i};
                $val =~s/\"//g;
                if ($val =~/\w+/){
                        $concat_accessory .= "\t1";
			$hash{$sample}{$cl_num} = 1;
                }
                else{
                        $concat_accessory .= "\t0";
			$hash{$sample}{$cl_num} = 0;
                }
        }
        if ($concat_accessory =~/0/){
                print M $cl_num.$concat_accessory."\n";
        }
	#if ($cl_num > 1000){last;}
}
close(F);
close(M);


system("ulimit -s 163840;Rscript $PANEX_PATH/R/heatmap.R -f $heatmap.accessory_01matrix.txt -o $heatmap.complete.svg");


my $min_y = 10000000;
my $min_x = 10000000;
my $max_y = 0;
my $max_x = 0;
my $step_x;
my $step_y;
open(SVG,"$heatmap.complete.svg");
while(<SVG>){
	if (/\<rect.*stroke/ && /x='(\d+\.\d+)' y='(\d+\.\d+)' width='(\d+\.\d+)' height='(\d+\.\d+)'/ ){
		my $y = $2;
		my $x = $1;
		$step_x = $3;
		$step_y = $4;
		#$y -= $height;
		if ($x < $min_x){$min_x = $x;}
		if ($y < $min_y){$min_y = $y;}
		if ($x > $max_x){$max_x = $x;}
		if ($y > $max_y){$max_y = $y;}
	}
}
close(SVG);


my @clusters;
open(F,"$heatmap.complete.svg.cols.csv");
<F>;
while(<F>){
	my $line = $_;
	$line =~s/\n//g;$line =~s/\r//g;
	push(@clusters,$line);
}
close(F);

my @strains;
open(F,"$heatmap.complete.svg.rows.csv");
<F>;
while(<F>){
        my $line = $_;
        $line =~s/\n//g;$line =~s/\r//g;
        push(@strains,$line);
}
close(F);

my $svg_section = "";

my %combinaisons;
my $x = $min_x;
my $length_matrix = ($max_x-$min_x);
my $height_matrix = ($max_y-$min_y);
my $step_x = $length_matrix / scalar @clusters;
#my $step_y = $height_matrix / scalar @strains;

print "$height_matrix $length_matrix $step_x $step_y\n";
my $previous_combinaison = "";
my $combinaison_id = 0;
my %combinaison_ids;
open(M,">$heatmap.upsetr.txt");
print M "ClutserID";
foreach my $strain(@strains){
	my @words = split(/_/,$strain);
        my $genus = $words[0];
        my $species = $words[1];
        my $shortname = substr($genus,0,3) . "_". substr($species,0,2);
        for (my $j = 2; $j <= $#words; $j++){
                $shortname.="_".$words[$j];
        }
        $shortname = substr($shortname,0,25);
	print M "\t".$shortname;
}
print M "\n";
foreach my $cluster(@clusters){
	$x = $x+$step_x;
	my $combinaison = "";
	my $numstrain;
	print M "$cluster";
	foreach my $strain(@strains){
		$numstrain++;
		my $val = $hash{$strain}{$cluster};
		print M "\t".$val;
		if ($val == 1){
			$combinaison .= ",$numstrain";
		}
	}
	if ($combinaison ne $previous_combinaison){
		$combinaison_id++;
	}
	print M "\n";
	$previous_combinaison = $combinaison; 
	$combinaison_ids{$combinaison_id}=$combinaison;
	#print "Cluster:$cluster Combinaison: $combinaison X:$x $combinaison_id\n";
	$combinaisons{$combinaison_id} .= "|".$x;
	
}
close(M);

my $nb_strains = scalar @strains;
system("Rscript $PANEX_PATH/R/upsetr.R $heatmap.upsetr.txt $heatmap.upsetr.pdf $nb_strains");
system("pdf2svg $heatmap.upsetr.pdf $heatmap.upsetr.svg 2");

my $numstrain;
my $y = $min_y + $height_matrix + $step_y;
foreach my $strain(@strains){
	$numstrain++;
	$y = $y - $step_y;
	my $x = $min_x;
	foreach my $combinaison_id(keys(%combinaisons)){
		my $combinaison = $combinaison_ids{$combinaison_id};
		my $x_values = $combinaisons{$combinaison_id};
		my @table_x = split(/\|/,$x_values);
		my $first_x = $table_x[1];
		my $last_x = $table_x[$#table_x];
		my $step_x = $last_x - $first_x;
		my $new_combinaison = "$combinaison,";
		if ($new_combinaison =~/,$numstrain,/){
			#$svg_section .= "<rect x=\"$first_x\" y=\"$y\" width=\"$step_x\" height=\"$step_y\" style=\"stroke: red; fill: red;\"><title>$strain $combinaison</title></rect>";
			$svg_section .= "<rect x=\"$first_x\" y=\"$y\" width=\"$step_x\" height=\"$step_y\" style=\"stroke: purple; fill: purple;\"/>";
		}
		#print "$combinaison $x_values $first_x $last_x\n";
	}
	#last;
}

open(FINAL_SVG,">$heatmap");
open(SVG,"$heatmap.complete.svg");
while(<SVG>){
	#if (/<path/){next;}
	if (/^\<rect/){next;}
	if (/^\<\/svg/){
		print FINAL_SVG $svg_section;
	}
	print FINAL_SVG $_;
}
close(SVG);
close(FINAL_SVG);
