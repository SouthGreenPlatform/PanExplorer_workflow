#!/usr/bin/perl
use strict;
use warnings;


my $BLASTResult=$ARGV[0];
my $SEQFile=$ARGV[1];
my $coverage=$ARGV[2];
my $identity=$ARGV[3]*100;#percents
my $score=$ARGV[4];
my %hash;
my @row;
my $line;

#my $coverage=0.5;
#my $identity=50;#percents
#my $score=50;

&ReadSeqLength("$SEQFile",\%hash);

open(F,$BLASTResult);
while ($line=<F>) 
	{
		if ($line!~/\#/) 
		{
			@row=split(/\t/,$line);
			if ($row[2]>=$identity and $row[11]>=$score) 
				{
					if ((($row[7]-$row[6]+1)/$hash{$row[0]})>=$coverage) 
						{
							if ((($row[9]-$row[8]+1)/$hash{$row[1]})>=$coverage) 
								{
									print "$row[0]\t$row[1]\t$row[11]";
								}
						}
				}
		}
	}
close(F);

sub ReadSeqLength()
	{
		use Bio::SeqIO;
		(my $file,my $hash)=@_;
		my $seq;
		my $in=Bio::SeqIO->new(-file=>"$file",-format=>"fasta");
		while ($seq=$in->next_seq()) 
			{
				$$hash{$seq->id}=length($seq->seq());
			}
	}