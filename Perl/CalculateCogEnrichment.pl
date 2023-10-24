#!/usr/bin/perl

use strict;

my $pav_matrix = $ARGV[0];
my $cog_of_clusters = $ARGV[1];
my $output = $ARGV[2];

my %categories = ("Unknown"=>1); 
my %cog_categories;
my %cogs_for_cluster;
open(COG,$cog_of_clusters);
while(<COG>){
	my $line = $_;
	$line =~s/\n//g;$line =~s/\r//g;
	my ($cluster,$cog,$category) = split("\t",$line);
	#if (!$cogs_for_cluster{$cluster}){
		$cogs_for_cluster{$cluster} .= "$cog,";
	#}
	$cog_categories{$cog} = $category;
	$categories{$category}=1;
}
close(COG);

my %hash;
my %hash2;
open(F,$pav_matrix);
my $first_line = <F>;
$first_line =~s/\n//g;$first_line =~s/\r//g;
my @samples = split(/\t/,$first_line);
my $nb_total = 0;
my $nb_core = 0;
my $nb_accessory = 0;
print C "Cluster\tFunction\tCOG\tCOG categories\n";
print D "Cluster\tFound in N strains\tFunction\tCOG\tCOG categories\n";
while(<F>){
	my @i = split("\t",$_);
	my $cluster_num = $i[0];
	my $clnb = $i[0];
	my $pos = $clnb;
	my $pos_before = $pos - 1;
	my $nb_found = 0;
	my $concatenate_samples = "";
	my $strain_found = 0;
	my %functions_of_genes;
	for (my $j = 1; $j <= $#i; $j++){
		my $val = $i[$j];
		if ($val =~/\w+/){
			my @genes = split(",",$val);
			$nb_found++;
		}
	}
	# core-genes
	my $type = "";
	if ($nb_found == (scalar @samples - 1)){
		$type = "core";
		$nb_core++;
	}
	# dispensables
	else{
		$type = "accessory";
		$nb_accessory++;
	}
	my @cogs = split(/,/,$cogs_for_cluster{$clnb});
	my %cat;
	foreach my $cog(@cogs){
		if ($cog =~/\w+/){
			my $category = $cog_categories{$cog};
			$cat{$category}=1;
		}
	}
	if (!$cogs_for_cluster{$clnb}){
		$cat{"Unknown"}=1;
	}
	
	foreach my $category(keys(%categories)){
		if ($cat{$category}){
			$hash{$category}{$type}++;
		}
		else{
			$hash{"not in $category"}{$type}++;
		}
	}
}
close(F);

open(OUT,">$output");
my %results;
foreach my $category(keys(%categories)){
	my $n1 = 0;
	if ($hash{$category}{"core"}){
		$n1 = $hash{$category}{"core"};
	}
	my $n2 = 0;
	if ($hash{$category}{"accessory"}){
		$n2 = $hash{$category}{"accessory"};
	}
	my $n3 = 0;
	if ($hash{"not in $category"}{"core"}){
		$n3 = $hash{"not in $category"}{"core"};
	}
	my $n4 = 0;
	if ($hash{"not in $category"}{"accessory"}){
		$n4 = $hash{"not in $category"}{"accessory"};
	}

	# core
	open(R,">$pav_matrix.$category.core.script.R");
	print R "(tab <- matrix(c($n1, $n2, $n3, $n4), nrow=2))\n";
	my $ratio1 = $n1/$n3;
	my $ratio2 = $n2/$n4;
	print OUT "core	$category	$n1	$n2	$n3	$n4	$ratio1	$ratio2\n";
	print R "fisher.test(tab)\n";
	close(R);
	system("Rscript $pav_matrix.$category.core.script.R >$pav_matrix.$category.core.script.R.out");
	open(O,"$pav_matrix.$category.core.script.R.out");
	my $go_odds = 0;
	while(<O>){
		if (/p-value [=<] (.*)/){
			my $pvalue = $1;
			$results{$category}{"core"}{"pval"} = $pvalue;
		}
		if (/odds/){$go_odds = 1;}
		if (/([\d\.]*)/ && $go_odds){
			my $odds = $1;
			
			$results{$category}{"core"}{"odds_ratio"} = $odds;
		}
	}
	close(O);

	if ($category ne 'Unknown'){
	unlink("$pav_matrix.$category.core.script.R.out");
	unlink("$pav_matrix.$category.core.script.R");
	}

	# dispensable
	my $n1 = 0;
	if ($hash{$category}{"accessory"}){
                $n1 = $hash{$category}{"accessory"};
        }
	my $n2 = 0;
	if ($hash{$category}{"core"}){
		$n2= $hash{$category}{"core"};
	}
	my $n3 = 0;
        if ($hash{"not in $category"}{"accessory"}){
                $n3 = $hash{"not in $category"}{"accessory"};
        }
	my $n4 = 0;
	if ($hash{"not in $category"}{"core"}){
		$n4 = $hash{"not in $category"}{"core"};
	}

	open(R,">$pav_matrix.$category.accessory.script.R");
        print R "(tab <- matrix(c($n1, $n3, $n2, $n4), nrow=2))\n";
        print R "fisher.test(tab)\n";
        close(R);
	print "dispensable, $category: $n1, $n2, $n3, $n4\n";
        system("Rscript $pav_matrix.$category.accessory.script.R >$pav_matrix.$category.accessory.script.R.out");
        open(O,"$pav_matrix.$category.accessory.script.R.out");
	my $go_odds = 0;
        while(<O>){
                if (/p-value [=<] (.*)/){
                        my $pvalue = $1;
                        $results{$category}{"accessory"}{"pval"} = $pvalue;
                }
		if (/odds/){$go_odds = 1;}
                if (/([\d\.]*)/ && $go_odds){
                        my $odds = $1;
                        $results{$category}{"accessory"}{"odds_ratio"} = $odds;
                }
        }
        close(O);
	if ($category ne 'Unknown'){
	unlink("$pav_matrix.$category.accessory.script.R");
	unlink("$pav_matrix.$category.accessory.script.R.out");
	}
}

foreach my $category(keys(%results)){
	print "$category:";
	print " ".$results{$category}{"core"}{"odds_ratio"};
	print " ".$results{$category}{"core"}{"pval"};
	print " ".$results{$category}{"accessory"}{"odds_ratio"};
	print " ".$results{$category}{"accessory"}{"pval"};
	print "\n";
}
close(OUT);
