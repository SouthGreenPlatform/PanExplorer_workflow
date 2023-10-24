#!/usr/bin/perl

use strict;

my $filein = $ARGV[0];
my $fileout = $ARGV[1];
my $matrix = $ARGV[2];




my $min_x = 100000;
my $max_x = 0;
my $min_y = 100000;
my $max_y = 0;
my $n=0;
open(OUT,">$fileout");
open(F,$filein);
while(<F>){
	if (/<path style=" stroke:none;fill-rule:nonzero;fill/){

		if (/d=\"M ([\s\d\.]+) L/){
			my @infos = split(/ /,$1);
			my $x = $infos[0];
			if ($x < $min_x){$min_x = $x;}
			if ($x > $max_x){$max_x = $x;}
		}
		if (/L ([\s\d\.]+) L/){
			my @infos = split(/ /,$1);
			my $y = $infos[1];
			if ($y < $min_y){$min_y = $y;}
			if ($y > $max_y){$max_y = $y;}
                }
		if (/L ([\s\d\.]+) Z/){
                        my @infos = split(/ /,$1);
                        my $y = $infos[1];
                        if ($y < $min_y){$min_y = $y;}
                        if ($y > $max_y){$max_y = $y;}
                }
		$n++;
	}
	else{
		if (!/\<\/svg\>/){
			print OUT $_;
		}
	}
}
close(F);

my $nb_dispensable_clusters = `grep -P -c '\t0' $matrix`;
my $nb_samples = `awk {'print NF-1'} $matrix | head -1`;

my $global_width = $max_x - $min_x;
my $width_of_one_block = $global_width / $nb_dispensable_clusters;

my $global_height = $max_y - $min_y;
my $height_of_one_block = $global_height / $nb_samples;

###########################################################
# get distinct pattern of presence/absence
###########################################################
my %patterns;
my $pattern_order = 0;
my %pattern_orders;
open(M,$matrix);
<M>;
while(<M>){
	my $line = $_;
	$line =~s/\n//g;$line =~s/\r//g;
	my @infos = split(/\t/,$line);
	
	my $pattern = "";
	for (my $k=1;$k<=$#infos;$k++){
		$pattern.=$infos[$k];
	}

	# print only dispensable (at least one absence)
	if ($pattern =~/0/){
		if (!$patterns{$pattern}){
			$pattern_order++;
		}
		$patterns{$pattern}++;
		$pattern_orders{$pattern_order} = $pattern;
	}
}
close(M);

print "Number of distinct patterns:";
print scalar keys(%patterns)."\n";

my @colors = ("orange","green","red","blue","black","pink","yellow","brown","grey","purple","darkred");

my $cumul_x = 0;
foreach my $pattern_order(sort {$a<=>$b} keys(%pattern_orders)){
	my $pattern = $pattern_orders{$pattern_order};
	my $size = $patterns{$pattern};
	my $width = $size * $width_of_one_block;
	my $x = $max_x - $cumul_x - $width;

	my $modulo = $pattern_order % 2;
	print "$pattern_order $pattern $size $modulo\n";
	
	#my $color = $colors[$pattern_order-1];
	my $color = $colors[$modulo];

	my $pattern_y = $min_y-15;
	print OUT "<rect y='$pattern_y' x='$x' width='$width' height='10' style=\"fill:$color;stroke-width:3;$color;\"/>";

	$cumul_x += $width;
	my @values = split(//,$pattern);
	my $cumul_y = 0;
	foreach my $val(@values){
		my $y = $max_y - $cumul_y - $height_of_one_block;
		if ($val){
			print OUT "<rect y='$y' x='$x' width='$width' height='$height_of_one_block' style=\"fill:purple;stroke:purple;\"/>";
		}
		$cumul_y += $height_of_one_block;
	}
}


print OUT "</svg>\n";
close(OUT);

print "Min x : $min_x\n";
print "Max x : $max_x\n";
print "Min y : $min_y\n";
print "Max y : $max_y\n";
print "$n\n";
