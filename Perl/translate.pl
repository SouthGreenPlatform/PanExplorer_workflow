#!/usr/bin/perl

use strict;

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
		system("perl /scratch2/galaxy/galaxy-20.09/galaxy/tools/SouthGreen/PanExplorer/PanExplorer_workflow/Perl/DNA-Transcription-Translation/DNA_Transcription_Translation.pl $DNA >>translate.log 2>&1");
		my $result = `cat result.txt`;
		chop($result);
		print O $result."\n";
	}
}
close(F);
close(O);
