#!/usr/bin/env python3
"""Keep one transcript per gene (the one with the longest CDS) in a GFF3 file.

Eukaryotic annotations often contain several transcripts (isoforms) per gene. In the PanExplorer
workflows every transcript would otherwise be treated as a separate "gene":
  - PGGB workflow: overlapping isoforms of the same gene compete for the same graph positions and are
    never linked together (genes of the same genome are not compared), which creates artificial
    singletons and inflates the accessory genome;
  - OrthoFinder workflow: all isoforms are clustered, which inflates gene counts per orthogroup.
Keeping the longest isoform in both workflows gives one entry per gene, so that the PAV matrices
of the two approaches can be compared gene to gene.

Transcripts are the features that are parents of CDS lines and whose own Parent is a gene. Lines of the
other transcripts of the same gene (mRNA, exon, CDS, UTR...) are removed. Genes without transcripts,
or CDS directly attached to a gene, are left untouched. A ##FASTA section is copied as is.

Usage: longest_isoform_gff.py input.gff3 output.gff3
"""
import sys
from collections import OrderedDict, defaultdict


def parse_attributes(col9):
    attrs = {}
    for field in col9.strip().split(";"):
        if "=" in field:
            key, value = field.split("=", 1)
            attrs.setdefault(key.strip(), value.strip())
    return attrs


def main(infile, outfile):
    header, records, fasta = [], [], []
    with open(infile) as handle:
        in_fasta = False
        for line in handle:
            if in_fasta:
                fasta.append(line)
                continue
            if line.startswith("##FASTA"):
                in_fasta = True
                fasta.append(line)
                continue
            if line.startswith("#") or not line.strip():
                (header if not records else records).append((line, None))
                continue
            cols = line.rstrip("\n\r").split("\t")
            if len(cols) < 9:
                records.append((line, None))
                continue
            records.append((line, cols))

    parent_of = {}                       # feature ID -> list of parent IDs
    cds_length = defaultdict(int)        # transcript ID -> CDS length
    for _, cols in records:
        if cols is None:
            continue
        attrs = parse_attributes(cols[8])
        parents = attrs.get("Parent", "").split(",") if attrs.get("Parent") else []
        if "ID" in attrs:
            parent_of[attrs["ID"]] = parents
        if cols[2] == "CDS":
            for parent in parents:
                cds_length[parent] += int(cols[4]) - int(cols[3]) + 1

    # transcripts grouped by gene (transcript = parent of CDS lines, having itself a parent)
    transcripts_of_gene = OrderedDict()
    for transcript in cds_length:
        genes = parent_of.get(transcript, [])
        if not genes:
            continue  # CDS attached directly to a gene, or transcript without gene: kept
        transcripts_of_gene.setdefault(genes[0], []).append(transcript)

    dropped = set()
    for gene, transcripts in transcripts_of_gene.items():
        if len(transcripts) < 2:
            continue
        best = max(transcripts, key=lambda t: cds_length[t])  # first one kept on ties
        dropped.update(t for t in transcripts if t != best)

    # remove dropped transcripts and all their descendants (exons, CDS, UTRs...)
    def is_dropped(feature_id, parents):
        if feature_id in dropped:
            return True
        stack, seen = list(parents), set()
        while stack:
            p = stack.pop()
            if p in dropped:
                return True
            if p in seen:
                continue
            seen.add(p)
            stack.extend(parent_of.get(p, []))
        return False

    kept_lines = removed_lines = 0
    with open(outfile, "w") as out:
        out.writelines(line for line, _ in header)
        for line, cols in records:
            if cols is not None:
                attrs = parse_attributes(cols[8])
                parents = attrs.get("Parent", "").split(",") if attrs.get("Parent") else []
                if is_dropped(attrs.get("ID"), parents):
                    removed_lines += 1
                    continue
                kept_lines += 1
            out.write(line)
        out.writelines(fasta)

    n_multi = sum(1 for t in transcripts_of_gene.values() if len(t) > 1)
    sys.stderr.write(
        f"[longest_isoform_gff] {infile}: {len(transcripts_of_gene)} genes with coding transcripts, "
        f"{n_multi} with several isoforms, {len(dropped)} isoforms removed "
        f"({removed_lines} lines removed, {kept_lines} kept)\n")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    main(sys.argv[1], sys.argv[2])
