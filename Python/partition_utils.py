#!/usr/bin/env python3
"""Helpers for the chromosome-partitioned PGGB workflow (Snakefile_pggb_heatmap_upset_COG).

  partition_utils.py fasta2gfa  group.fa group.gfa
      GFA for a group containing a single genome (nothing to align): one segment and one path per sequence.
  partition_utils.py merge-pav  out.tsv groups_out.tsv group1:pav1.tsv group2:pav2.tsv ...
      concatenates the PAV matrices of the groups (same columns) and renumbers the clusters;
      groups_out.tsv gives the group (chromosome) of each cluster.
  partition_utils.py merge-vcf  out.vcf in1.vcf in2.vcf ...
      concatenates VCF files (header of the first non-empty file).
  partition_utils.py png  out.png
      writes a 1x1 white PNG (placeholder when no group produced a PGGB image).
"""
import base64
import sys


def fasta2gfa(fasta, gfa):
    with open(fasta) as handle, open(gfa, "w") as out:
        out.write("H\tVN:Z:1.0\n")
        names, seq, node = [], [], 0

        def flush():
            nonlocal node
            if names and seq:
                node += 1
                out.write(f"S\t{node}\t{''.join(seq)}\n")
                out.write(f"P\t{names[-1]}\t{node}+\t*\n")

        for line in handle:
            if line.startswith(">"):
                flush()
                names.append(line[1:].split()[0])
                seq = []
            else:
                seq.append(line.strip())
        flush()


def merge_pav(out, groups_out, items):
    header, n = None, 0
    with open(out, "w") as fout, open(groups_out, "w") as gout:
        gout.write("ClutserID\tgroup\n")
        for item in items:
            group, path = item.split(":", 1)
            with open(path) as handle:
                head = handle.readline()
                if not head:
                    continue
                if header is None:
                    header = head
                    fout.write(head)
                elif head != header:
                    sys.exit(f"[merge-pav] columns of {path} differ from the first matrix")
                for line in handle:
                    if not line.strip():
                        continue
                    n += 1
                    fout.write(str(n) + line[line.index("\t"):])
                    gout.write(f"{n}\t{group}\n")
    sys.stderr.write(f"[merge-pav] {n} clusters from {len(items)} groups\n")


def merge_vcf(out, inputs):
    header_done = False
    with open(out, "w") as fout:
        for path in inputs:
            with open(path) as handle:
                lines = handle.readlines()
            if not lines:
                continue
            for line in lines:
                if line.startswith("#"):
                    if not header_done:
                        fout.write(line)
                else:
                    fout.write(line)
            header_done = True


def png(out):
    data = ("iVBORw0KGgoAAAANSUhEUgAAAAEAAAABCAYAAAAfFcSJAAAADUlEQVR42mP8/58BAwAI/AL+hc2rNAAAAABJRU5ErkJggg==")
    with open(out, "wb") as fout:
        fout.write(base64.b64decode(data))


if __name__ == "__main__":
    if len(sys.argv) < 3:
        sys.exit(__doc__)
    cmd = sys.argv[1]
    if cmd == "fasta2gfa":
        fasta2gfa(sys.argv[2], sys.argv[3])
    elif cmd == "merge-pav":
        merge_pav(sys.argv[2], sys.argv[3], sys.argv[4:])
    elif cmd == "merge-vcf":
        merge_vcf(sys.argv[2], sys.argv[3:])
    elif cmd == "png":
        png(sys.argv[2])
    else:
        sys.exit(__doc__)
