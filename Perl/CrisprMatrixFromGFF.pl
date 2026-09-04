#!/usr/bin/env perl

use strict;
use warnings;
use File::Basename;

# ------------------------------------------------------------
# Usage:
#   perl gff_crispr_matrix.pl /path/to/gff_directory output.tsv
#
# Example:
#   perl gff_crispr_matrix.pl ./gff crispr_presence_absence.tsv
# ------------------------------------------------------------

my $gff_dir = shift @ARGV
    or die "Usage: $0 <gff_directory> <output.tsv>\n";

my $output = shift @ARGV
    or die "Usage: $0 <gff_directory> <output.tsv>\n";

# ------------------------------------------------------------
# Find GFF files
# ------------------------------------------------------------

opendir(my $DIR, $gff_dir)
    or die "Cannot open directory $gff_dir: $!\n";

my @files = grep {
       /\.gff3?$/i
    || /\.gff3?\.gz$/i
} readdir($DIR);

closedir($DIR);

die "No GFF files found in $gff_dir\n"
    unless @files;

@files = sort @files;

# ------------------------------------------------------------
# Data structures
#
# %presence
#   $presence{spacer_sequence}{genome} = 1
#
# %spacer_id
#   $spacer_id{sequence} = SP001, SP002, ...
# ------------------------------------------------------------

my %presence;
my %spacer_id;

my $spacer_number = 0;

# ------------------------------------------------------------
# Read each GFF
# ------------------------------------------------------------

foreach my $file (@files) {

    my $path = "$gff_dir/$file";

    # Genome name = filename without extensions
    my $genome = basename($file);
    $genome =~ s/\.gz$//i;
    $genome =~ s/\.gff3?$//i;

    print STDERR "Reading $file ...\n";

    my $FH;

    if ($file =~ /\.gz$/i) {
        open($FH, "gzip -dc '$path' |")
            or die "Cannot open $path: $!\n";
    }
    else {
        open($FH, "<", $path)
            or die "Cannot open $path: $!\n";
    }

    while (my $line = <$FH>) {

        chomp $line;

        # Skip comments
        next if $line =~ /^#/;

        my @fields = split(/\t/, $line);

        next unless @fields >= 9;

        my $type = $fields[2];

        # We only want CRISPR spacers
        next unless $type eq "CRISPRspacer";

        my $attributes = $fields[8];

        # Extract sequence=...
        if ($attributes =~ /(?:^|;)sequence=([^;]+)/) {

            my $sequence = uc($1);

            # Remove possible whitespace
            $sequence =~ s/\s+//g;

            # Basic DNA validation
            unless ($sequence =~ /^[ACGTN]+$/) {
                warn "Warning: unusual spacer sequence in $file: $sequence\n";
            }

            # Assign a global spacer ID
            unless (exists $spacer_id{$sequence}) {
                $spacer_number++;

                my $id = sprintf("SP%04d", $spacer_number);

                $spacer_id{$sequence} = $id;
            }

            # Presence = 1
            $presence{$sequence}{$genome} = 1;
        }
    }

    close($FH);
}

# ------------------------------------------------------------
# Sort spacers by their assigned ID
# ------------------------------------------------------------

my @spacers = sort {
    $spacer_id{$a} cmp $spacer_id{$b}
} keys %spacer_id;

my @genomes;

foreach my $file (@files) {

    my $genome = basename($file);
    $genome =~ s/\.gz$//i;
    $genome =~ s/\.gff3?$//i;

    push @genomes, $genome;
}

# ------------------------------------------------------------
# Write matrix
# ------------------------------------------------------------

open(my $OUT, ">", $output)
    or die "Cannot write $output: $!\n";

# Header
print $OUT "Spacer_ID\tSequence\t", join("\t", @genomes), "\n";

# Rows
foreach my $sequence (@spacers) {

    my $id = $spacer_id{$sequence};

    print $OUT $id, "\t", $sequence;

    foreach my $genome (@genomes) {

        my $value = exists $presence{$sequence}{$genome}
                  ? 1
                  : 0;

        print $OUT "\t", $value;
    }

    print $OUT "\n";
}

close($OUT);

# ------------------------------------------------------------
# Summary
# ------------------------------------------------------------

print STDERR "\n";
print STDERR "========================================\n";
print STDERR "CRISPR spacer analysis\n";
print STDERR "========================================\n";
print STDERR "GFF files      : ", scalar(@files), "\n";
print STDERR "Genomes        : ", scalar(@genomes), "\n";
print STDERR "Unique spacers : ", scalar(@spacers), "\n";
print STDERR "Output         : $output\n";
print STDERR "========================================\n";
