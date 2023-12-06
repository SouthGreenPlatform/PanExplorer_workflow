#!/usr/bin/perl

use strict;

my $pav = $ARGV[0];
my $ref = $ARGV[1];
my $strain_info = $ARGV[2];
my $gff = $ARGV[3];
my $ordering = $ARGV[4];
my $outfile = $ARGV[5];

my $circos_exe = "circos/bin/circos";
#my $circos_exe = "circos";

my %metadata;
open(M,$strain_info);
while(my $line = <M>){
	chomp($line);
	my ($id,$name) = split(/\t/,$line);
	$metadata{$id} = $name;
	$metadata{$name} = $id;
}
close(M);


my $pangenome_size = 0;
my $nb_core_genes = 0;
my $nb_dispensable_genes = 0;
my $nb_specific_genes = 0;
my %cluster_of_gene;
my %genes_of_cluster;
my %tagging;
open(P,$pav);
my $firstline = <P>;
chomp($firstline);
my @infos = split(/\t/,$firstline);
my %samples;
my %hash_presence;
my $col;
for (my $j = 1; $j <= $#infos; $j++){
	my $val = $infos[$j];
	if ($val eq $ref){$col = $j;}
	$samples{$j} = $val;
}
while(my $line = <P>){
	chomp($line);
	my @infos = split(/\t/,$line);
	my $cluster_num = $infos[0];
        my $clnb = $infos[0];
        my $pos = $clnb;
        my $pos_before = $pos - 1;
        $pangenome_size++;
        if (length($clnb) == 1){$clnb = "000".$clnb;}
        elsif (length($clnb) == 2){$clnb = "00".$clnb;}
        elsif (length($clnb) == 3){$clnb = "0".$clnb;}
        my $name = "CLUSTER".$clnb;

	# exclude lines where several genes for one strains
	if ($line =~/,/){next;}

	my $gene_ref = $infos[$col];
	for (my $j = 1; $j <= $#infos; $j++){
		my $gene = $infos[$j];
		my $sample = $samples{$j};
		if ($gene =~/\w+/){
			$hash_presence{$gene_ref}{$sample} = 1;
		}
	}
}
close(P);

my $track_thin = 0.5 / scalar keys(%samples);
my $r0 = 0.4;
my $r1 = $r0+$track_thin;
my $highlight_block = "";
my %max_position;
open(O,$ordering);
<O>;
while(my $sample_name = <O>){
#foreach my $sample(keys(%samples)){
	chomp($sample_name);
	#my $sample_name = $samples{$sample};
	my $sample_id = $metadata{$sample_name};
	print "$sample_name $sample_id\n";
	my $ref_id = $metadata{$ref};
	open(GFF,$gff);
	open(CIRCOS_TRACK,">$outfile.$sample_id.circos.heatmap.txt");
	my $current_start;
	my $current_stop;
	my $current_chrom;
	$r0+=$track_thin;
	$r1+=$track_thin;
	#my $r1_final = $r1-($track_thin/3);
	my $r1_final = $r1;
	$r0.="r";
	$r1_final.="r";	
	$highlight_block .= qq~
<highlight>
file       = $outfile.$sample_id.circos.heatmap.txt
r0         = $r0
r1         = $r1_final
</highlight>
~;
	while(my $line = <GFF>){
		chomp($line);
		my @infos = split(/\t/,$line);
		if ($infos[2] eq "gene" && $line =~/ID=([^;]+);*/){
			my $chr = $infos[0];
			my $start = $infos[3];
			my $end = $infos[4];
			$current_start = $start;
			$current_stop = $end;
			$current_chrom = $chr;
			if ($end > $max_position{$current_chrom}){$max_position{$current_chrom} = $end;}
			my $gene = $1;
			my $genelength = $end-$start;
			my $gene_complete = $sample_name.":".$gene;
			if ($hash_presence{$gene}{$sample_name}){
			}	
		}
		if ($infos[2] eq "CDS" && $line =~/Parent=([^;]+);*/){
			my $gene = $1;
			my $gene_complete = $sample_name.":".$gene;
			if ($line =~/protein_id=([^;]+);/){
				$gene = $1;
				if ($hash_presence{$gene}{$sample_name}){
					print CIRCOS_TRACK "$current_chrom $current_start $current_stop fill_color=red\n";
				}
			}
		}	
	}
	close(GFF);
	close(CIRCOS_TRACK);
}
open(KARYOTYPE,">$outfile.karyotype.txt");
foreach my $chr(keys(%max_position)){
	my $max = $max_position{$chr};
	print KARYOTYPE "chr - $chr 1 0 $max black";
}
close(KARYOTYPE);

open(F,"circos_templates/circos1.conf");
open(CIRCOS_CONF,">$outfile.circos1.conf");
while(<F>){
	if (/HIGHLIGHT_BLOCK/){print CIRCOS_CONF $highlight_block;}
	elsif (/KARYOTYPE_PATH/){print CIRCOS_CONF "karyotype = $outfile.karyotype.txt\n";}
	elsif (/OUTPUT_PNG/){print CIRCOS_CONF "file  = $outfile\n";}
	else{print CIRCOS_CONF $_;}
}
close(F);
close(CIRCOS_CONF);

system("$circos_exe -conf $outfile.circos1.conf");
