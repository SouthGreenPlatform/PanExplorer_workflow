#!/usr/bin/perl

use strict;
use File::Basename;
my $dirname = dirname(__FILE__);

my $infile = $ARGV[0];
my $heatmap = $ARGV[1];

my %core_cluster;
my %samples;
my %sequences;
my %matrix;
my %hash;
my $cl_num = 0;
my $nb_strains = 1;
open(M,">$heatmap.accessory_01matrix.txt");
open(PAN,">$heatmap.pangenes_01matrix.txt");
open(F,"$infile");
my $firstline = <F>;
$firstline =~s/\n//g;$firstline =~s/\r//g;
#my @infos = split(/","/,$firstline);
my @infos = split(/\t/,$firstline);
print M "Gene";
print PAN "Gene";
#for (my $j=14; $j <= $#infos; $j++){
for (my $j=1; $j <= $#infos; $j++){
        my $gbfile = $infos[$j];
        $gbfile =~s/\"//g;
        $gbfile =~s/\.gb\.filt//g;
	
	my $strain = $gbfile;
	my @words = split(/_/,$strain);
        my $genus = $words[0];
        my $species = $words[1];
        my $shortname = substr($genus,0,3) . "_". substr($species,0,2);
        for (my $j = 2; $j <= $#words; $j++){
                $shortname.="_".$words[$j];
        }
        #$shortname = substr($shortname,0,25);

        print M "\t".$strain;
	print PAN "\t".$strain;
	$samples{$j} = $strain;
        $nb_strains++;
}
print M "\n";
print PAN "\n";
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
		$matrix{$sample}{$cl_num} = $val;
                $val =~s/\"//g;
                if ($val =~/\w+/){
                        $concat_accessory .= "\t1";
			$hash{$sample}{$cl_num} = 1;
			$sequences{$sample}.=  "T";
                }
                else{
                        $concat_accessory .= "\t0";
			$hash{$sample}{$cl_num} = 0;
			$sequences{$sample}.=  "A";
                }
        }
        if ($concat_accessory =~/0/){
                print M $cl_num.$concat_accessory."\n";
        }
	else{
                $core_cluster{$cl_num}=1;
        }
	print PAN $cl_num.$concat_accessory."\n";
	#if ($cl_num > 1000){last;}
}
close(F);
close(M);
close(PAN);


system("ulimit -s 163840;Rscript $dirname/../R/heatmap.R -f $heatmap.accessory_01matrix.txt -o $heatmap.complete.pdf");
#system("ulimit -s 163840;Rscript $dirname/../R/heatmaply.R -f $heatmap.accessory_01matrix.txt -o $heatmap.complete2.svg");
system("pdf2svg $heatmap.complete.pdf $heatmap.complete.svg");


open(A,">$heatmap.alignment.fa");
foreach my $sample(keys(%sequences)){
	my $seq = $sequences{$sample};
	print A ">$sample\n";
	print A "$seq\n";
}
close(A);



my $min_y = 10000000;
my $min_x = 10000000;
my $max_y = 0;
my $max_x = 0;
my $step_x;
my $step_y;
open(SVG,"$heatmap.complete.svg");
open(SVGNEW,">$heatmap.complete.new.svg");
while(<SVG>){
	if (/rgb\(100%,0%,0%\)/){
		$_ =~s/rgb\(100%,0%,0%\)/rgb\(100%,100%,100%\)/g;
	}
	if (/rgb\(100%,63.529968%,0%\)/){
                $_ =~s/rgb\(100%,63.529968%,0%\)/rgb\(41%,72%,64%\)/g;
        }
	print SVGNEW $_;
}
close(SVG);
close(SVGNEW);



my @clusters;
open(F,"$heatmap.complete.pdf.cols.csv");
<F>;
while(<F>){
	my $line = $_;
	$line =~s/\n//g;$line =~s/\r//g;
	push(@clusters,$line);
	print $line;
}
close(F);

my @strains;
open(F,"$heatmap.complete.pdf.rows.csv");
<F>;
while(<F>){
        my $line = $_;
        $line =~s/\n//g;$line =~s/\r//g;
        push(@strains,$line);
	print $line;
}
close(F);


open(O,">$infile.sorted");
open(HP,">$infile.sorted.for_heatmap_plotly.txt");
print O "ClutserID"."\t".join("\t",@strains)."\n";
print HP "ClutserID"."\t".join("\t",@strains)."\n";
foreach my $cl(reverse(@clusters)){
	print O $cl;
	my $clnb = $cl;
	if (length($clnb) == 1){$clnb = "000".$clnb;}
        elsif (length($clnb) == 2){$clnb = "00".$clnb;}
        elsif (length($clnb) == 3){$clnb = "0".$clnb;}
        my $name = "CLUSTER".$clnb;

	print HP $name;
        foreach my $strain(@strains){
                my $val = $matrix{$strain}{$cl};
                print O "\t$val";
		if ($val =~/\w/){print HP "\t1";}
		else{print HP "\t0";}
        }
        print O "\n";
	print HP "\n";
}
foreach my $cl(keys(%core_cluster)){
        print O $cl;
	my $clnb = $cl;
	if (length($clnb) == 1){$clnb = "000".$clnb;}
        elsif (length($clnb) == 2){$clnb = "00".$clnb;}
        elsif (length($clnb) == 3){$clnb = "0".$clnb;}
        my $name = "CLUSTER".$clnb;
	print HP $name;
        foreach my $strain(@strains){
                my $val = $matrix{$strain}{$cl};
                print O "\t$val";
		if ($val =~/\w/){print HP "\t1";}
		else{print HP "\t0";}
		
        }
        print O "\n";
	print HP "\n";
}
close(O);
close(HP);

my $svg_section = "";

if (scalar @clusters == 0){
	system("touch $heatmap.upsetr.svg");
	if (!-e "$heatmap"){
		system("touch $heatmap");
	}	
	exit;
}
system("cp -rf $infile $infile.bkp");
system("cp -rf $infile.sorted $infile");

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
system("Rscript $dirname/../R/upsetr.R $heatmap.upsetr.txt $heatmap.upsetr.pdf $nb_strains >>$heatmap.upsetr.log 2>&1");
system("pdf2svg $heatmap.upsetr.pdf $heatmap.upsetr.svg 2");

system("python3 $dirname/../Python/Heatmap.py -i $infile.sorted.for_heatmap_plotly.txt -o $heatmap.heatmap_plotly.html >>$heatmap.heatmap_plotly.log 2>&1");

if (!-e "$heatmap.upsetr.svg"){
        system("touch $heatmap.upsetr.svg");
}
system("perl $dirname/../Perl/reformatHeatmapSVG.pl $heatmap.complete.svg $heatmap.complete.new.svg $infile.sorted.for_heatmap_plotly.txt");
system("cp -rf $heatmap.complete.new.svg $heatmap");
system("gzip $heatmap");
