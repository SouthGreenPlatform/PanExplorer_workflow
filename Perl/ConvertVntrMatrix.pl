#!/usr/bin/perl

use strict;

my $matrix = $ARGV[0];
my $out = $ARGV[1];
my $strain_names = $ARGV[2];

my %strains_of_gb;
open(F,$strain_names);
while(<F>){
	my $line = $_;
	$line =~s/\n//g;$line =~s/\r//g;
	my ($gb,$strain) = split(/\t/,$line);
	$strains_of_gb{$gb} = $strain;
}
close(F);


my $nb_strains = 1;
open(O,">$out");
open(F,$matrix);
my $first_line = <F>;
print O $first_line;
my $n = 0;
my $monomorphic = 0;
my $m = 0;
while(my $line = <F>){
	chomp($line);
	my @infos = split(/\t/,$line);
       	my %hash;
	my %alleles;
	my $concat = $infos[0];	
	for (my $i = 1; $i <= $#infos; $i++){
		if ($infos[$i] !~/-/){
			my ($repeat_name,$repeat,$nb_repeat) = split(/_/,$infos[$i]);
			if ($nb_repeat =~/(\d+)\.+/){$nb_repeat = $1;}
			my $size_repeat = length($repeat);
			$concat .= "\t".$repeat."($nb_repeat)";
			$alleles{$nb_repeat}++;
			$hash{$size_repeat}++;
		}
		else{
			$concat .= "\t".$infos[$i];
		}
	}
	# different patterns
	if (scalar keys(%hash) > 1){
		$m++;
	}
	# monomorphic
	elsif (scalar keys(%alleles) == 1){
		$monomorphic++;
	}
	# found several times 
	elsif ($line =~/,/){
		
	}
	else{
		print O "$concat\n";
		$n++;
	}
}
close(F);
close(O);

print "$n $m $monomorphic\n";
