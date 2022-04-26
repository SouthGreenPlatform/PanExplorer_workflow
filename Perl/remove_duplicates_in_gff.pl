#!/usr/bin/perl

use strict;

my $file = $ARGV[0];
my $out = $ARGV[1];

my %h;
open(F,$file);
while(<F>){
	my @infos = split(/\t/,$_);
	if ($infos[2] eq 'CDS' && /Name=([^;]*);/){
		my $id = $1;
		$h{$id}++;
	}
}
close(F);

my %dup;
foreach my $id(keys(%h)){
	my $n = $h{$id};
	if ($n > 1){
		$dup{$id} = 1;
	}
}
open(O,">$out");
open(F,$file);
while(<F>){
        my @infos = split(/\t/,$_);
        if ($infos[2] eq 'CDS' && /Name=([^;]*);/){
		my $id = $1;
		if ($dup{$id}){next;}
	}
	if ($infos[2] eq 'mRNA' && /Parent=([^;]*);/){
                my $id = $1;
                if ($dup{$id}){next;}
        }
	if ($infos[2] eq 'gene' && /ID=([^;]*);/){
                my $id = $1;
                if ($dup{$id}){next;}
        }
	if ($infos[2] eq 'exon' && /Parent=([^;]*);/){
                my ($id,$extension) = split(/\./,$1);
                if ($dup{$id}){next;}
        }
	print O $_;

}
close(F);
close(O);
