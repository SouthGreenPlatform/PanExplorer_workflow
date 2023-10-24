#!/usr/bin/perl

use strict;

use File::Basename;
my $dirname = dirname(__FILE__);

my $out = $ARGV[1];
my $in = $ARGV[0];
open(F,$in);
open(O,">$out");
while(<F>){
	if (/>(.*)/){
		print O  $_;
	}
	else{
		my $dna = $_;
		my $DNA = uc($dna);
		$DNA =~s/\n//g;$DNA =~s/\r//g;
		system("perl $dirname/DNA_Transcription_Translation.pl $DNA >>translate.log 2>&1");
		my $result = `cat result.txt`;
		chop($result);
		print O $result."\n";
	}
}
close(F);
close(O);
