package GeneId;

# Shared rules for naming genes consistently across the PanExplorer workflow
# (GFF parsing for the PGGB gene projection, gene lengths, FASTA headers, matrix conversion).
#
# Prokaryotic annotations (no mRNA features: NCBI GenBank/RefSeq, Prokka, Bakta...):
#     locus_tag  >  protein_id  >  ID  >  Parent
#   locus_tag is unique for each gene, whereas RefSeq protein_id (WP_...) is shared by all
#   identical proteins, including multi-copy genes of the same genome (IS, TAL effectors...).
#   CDS split over several lines (join(), programmed frameshifts) share the same locus_tag
#   and are therefore merged into a single gene.
#
# Eukaryotic annotations (mRNA features present): unchanged behaviour,
#     protein_id  >  Parent  >  ID
#   so that one entry per transcript/protein is kept, matching the transcript-based
#   sequence files used downstream.

use strict;
use warnings;
use Exporter 'import';
our @EXPORT_OK = qw(gff_is_eukaryotic gene_key_from_attributes parse_gff_attributes);

# Returns 1 if the GFF file contains mRNA features (eukaryotic gene models), 0 otherwise.
sub gff_is_eukaryotic {
    my ($file) = @_;
    open(my $fh, '<', $file) or return 0;
    while (my $line = <$fh>) {
        next if $line =~ /^#/;
        my @c = split(/\t/, $line);
        if (defined $c[2] && $c[2] eq 'mRNA') {
            close($fh);
            return 1;
        }
    }
    close($fh);
    return 0;
}

# Parses the 9th column of a GFF3 line into a hash (first value kept for each key).
# Keys are matched exactly, so "old_locus_tag" is never confused with "locus_tag".
sub parse_gff_attributes {
    my ($attributes) = @_;
    my %a;
    return \%a unless defined $attributes;
    $attributes =~ s/[\r\n]+$//;
    foreach my $kv (split(/;/, $attributes)) {
        $kv =~ s/^\s+|\s+$//g;
        my ($k, $v) = split(/=/, $kv, 2);
        next unless defined $k && defined $v && $v ne '';
        $a{$k} = $v unless exists $a{$k};
    }
    return \%a;
}

# Returns the gene identifier to use for a CDS line, given its attribute column.
sub gene_key_from_attributes {
    my ($attributes, $is_eukaryotic) = @_;
    my $a = parse_gff_attributes($attributes);
    my @order = $is_eukaryotic ? qw(protein_id Parent ID)
                               : qw(locus_tag protein_id ID Parent);
    foreach my $key (@order) {
        return $a->{$key} if defined $a->{$key};
    }
    return undef;
}

1;
