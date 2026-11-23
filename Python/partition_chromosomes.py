#!/usr/bin/env python3
"""Group the chromosomes of several genomes into homologous groups (chr1 with chr1, ...) before PGGB.

Input: a PanSN FASTA (headers "genome#hap#sequence", as built by the PanExplorer workflows).
Each sequence of the other genomes is mapped on the sequences of a reference genome with an approximate
wfmash mapping (no base-level alignment, fast) and assigned to the reference sequence on which it has the
largest mapped length. The groups therefore do not depend on sequence names ("1", "chr01", "Chr1"...).

  - reference sequences shorter than --min-ref-length (unplaced scaffolds, organelles...) do not make their
    own group: they go to the group "unplaced", with the sequences mapping on them;
  - a sequence whose mapped length is below --min-frac of its length goes to "unplaced";
  - if wfmash fails or maps nothing, the sequences are grouped by normalised name (chr01 = Chr1 = 1).

Outputs (in --outdir):
  sequences.tsv   one line per sequence: PanSN name, genome, sequence name, length, group, method,
                  best reference sequence, fraction of the sequence mapped on it
  groups.tsv      one line per group: group, number of genomes, number of sequences, total length,
                  genomes (comma separated), reference present (yes/no) - sorted by decreasing length
  fasta/<group>.fa  sequences of each group
  partition_report.txt  human-readable summary with warnings

Usage: partition_chromosomes.py --fasta all_genomes.fa --reference GENOME --outdir outputs/partition
                                [--threads 4] [--min-frac 0.2] [--min-ref-length 1000000] [--paf FILE]
"""
import argparse
import collections
import os
import re
import shutil
import subprocess
import sys

DELIM = "#"


def read_fasta_index(fasta):
    """Return the list of (name, length) in file order (from the sequences, no faidx needed)."""
    seqs, name, length = [], None, 0
    with open(fasta) as handle:
        for line in handle:
            if line.startswith(">"):
                if name is not None:
                    seqs.append((name, length))
                name, length = line[1:].split()[0], 0
            else:
                length += len(line.strip())
    if name is not None:
        seqs.append((name, length))
    return seqs


def split_pansn(name):
    parts = name.split(DELIM, 2)
    if len(parts) == 3:
        return parts[0], parts[2]
    return name, name


def normalise_name(seq):
    """chr01 / Chr1 / CHR_1 / 1 -> '1' (used only when the mapping gives nothing)."""
    s = re.sub(r"^(chromosome|chrom|chr)[_\-\s]*", "", seq, flags=re.I)
    s = re.sub(r"\.\d+$", "", s)
    s = s.lstrip("0") or s
    return s.upper()


def safe_group_name(name):
    return re.sub(r"[^A-Za-z0-9_.-]", "_", name)


def write_subset(fasta, names, out):
    keep = set(names)
    with open(fasta) as handle, open(out, "w") as fout:
        write = False
        for line in handle:
            if line.startswith(">"):
                write = line[1:].split()[0] in keep
            if write:
                fout.write(line)


def run_wfmash(ref_fa, query_fa, paf, threads, identity):
    if shutil.which("wfmash") is None:
        sys.stderr.write("[partition] wfmash not found\n")
        return False
    for fa in (ref_fa, query_fa):
        subprocess.run(["samtools", "faidx", fa], check=False)
    cmd = ["wfmash", "--approx-map", "-s", "10000", "-l", "50000", "-p", str(identity), "-n", "1",
           "-t", str(threads), ref_fa, query_fa]
    sys.stderr.write("[partition] " + " ".join(cmd) + "\n")
    with open(paf, "w") as out:
        proc = subprocess.run(cmd, stdout=out)
    return proc.returncode == 0 and os.path.getsize(paf) > 0


def parse_paf(paf):
    """query -> target -> mapped query length (bp)."""
    mapped = collections.defaultdict(lambda: collections.defaultdict(int))
    with open(paf) as handle:
        for line in handle:
            f = line.split("\t")
            if len(f) < 12:
                continue
            mapped[f[0]][f[5]] += int(f[3]) - int(f[2])
    return mapped


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--reference", required=True, help="genome (PanSN prefix) used to define the groups")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--min-frac", type=float, default=0.2)
    ap.add_argument("--min-ref-length", type=int, default=1000000)
    ap.add_argument("--identity", type=float, default=85)
    ap.add_argument("--paf", help="use this mapping (query = other genomes, target = reference) instead of running wfmash")
    args = ap.parse_args()

    os.makedirs(os.path.join(args.outdir, "fasta"), exist_ok=True)
    seqs = read_fasta_index(args.fasta)
    length = dict(seqs)
    genome_of = {n: split_pansn(n)[0] for n, _ in seqs}
    seqname_of = {n: split_pansn(n)[1] for n, _ in seqs}
    genomes = list(dict.fromkeys(genome_of[n] for n, _ in seqs))
    if args.reference not in genomes:
        sys.exit(f"[partition] reference genome '{args.reference}' not found in {args.fasta} (genomes: {', '.join(genomes)})")

    ref_seqs = [n for n, _ in seqs if genome_of[n] == args.reference]
    query_seqs = [n for n, _ in seqs if genome_of[n] != args.reference]

    # groups defined by the reference sequences
    group_of_ref = {}
    for n in ref_seqs:
        group_of_ref[n] = safe_group_name(seqname_of[n]) if length[n] >= args.min_ref_length else "unplaced"

    assignment = {}  # seq -> (group, method, best_target, fraction)
    for n in ref_seqs:
        assignment[n] = (group_of_ref[n], "reference", n, 1.0)

    mapped = {}
    paf = args.paf
    if paf is None and query_seqs:
        work = os.path.join(args.outdir, "mapping")
        os.makedirs(work, exist_ok=True)
        ref_fa, query_fa = os.path.join(work, "reference.fa"), os.path.join(work, "query.fa")
        write_subset(args.fasta, ref_seqs, ref_fa)
        write_subset(args.fasta, query_seqs, query_fa)
        paf = os.path.join(work, "query_vs_reference.paf")
        if not run_wfmash(ref_fa, query_fa, paf, args.threads, args.identity):
            paf = None
    if paf:
        mapped = parse_paf(paf)

    name_index = {}
    for n in ref_seqs:
        name_index.setdefault(normalise_name(seqname_of[n]), n)

    for n in query_seqs:
        hits = mapped.get(n, {})
        if hits:
            best, bp = max(hits.items(), key=lambda kv: kv[1])
            frac = bp / length[n] if length[n] else 0
            if frac >= args.min_frac and best in group_of_ref:
                assignment[n] = (group_of_ref[best], "mapping", best, round(frac, 3))
                continue
            assignment[n] = ("unplaced", "mapping_low", best, round(frac, 3))
            continue
        # no mapping at all: fall back on names
        ref = name_index.get(normalise_name(seqname_of[n]))
        if ref is not None:
            assignment[n] = (group_of_ref[ref], "name", ref, 0.0)
        else:
            assignment[n] = ("unplaced", "none", "-", 0.0)

    # write the sequences of each group
    by_group = collections.OrderedDict()
    for n, _ in seqs:
        by_group.setdefault(assignment[n][0], []).append(n)
    handles = {g: open(os.path.join(args.outdir, "fasta", g + ".fa"), "w") for g in by_group}
    with open(args.fasta) as handle:
        out = None
        for line in handle:
            if line.startswith(">"):
                out = handles[assignment[line[1:].split()[0]][0]]
            out.write(line)
    for h in handles.values():
        h.close()

    with open(os.path.join(args.outdir, "sequences.tsv"), "w") as out:
        out.write("pansn_name\tgenome\tsequence\tlength\tgroup\tmethod\tbest_reference_sequence\tmapped_fraction\n")
        for n, _ in seqs:
            g, method, best, frac = assignment[n]
            out.write(f"{n}\t{genome_of[n]}\t{seqname_of[n]}\t{length[n]}\t{g}\t{method}\t{best}\t{frac}\n")

    rows = []
    for g, members in by_group.items():
        gens = list(dict.fromkeys(genome_of[n] for n in members))
        rows.append((g, len(gens), len(members), sum(length[n] for n in members), ",".join(gens),
                     "yes" if args.reference in gens else "no"))
    rows.sort(key=lambda r: -r[3])
    with open(os.path.join(args.outdir, "groups.tsv"), "w") as out:
        out.write("group\tn_genomes\tn_sequences\ttotal_length\tgenomes\treference_present\n")
        for r in rows:
            out.write("\t".join(map(str, r)) + "\n")

    with open(os.path.join(args.outdir, "partition_report.txt"), "w") as rep:
        rep.write(f"Reference genome: {args.reference}\nGenomes: {', '.join(genomes)}\n")
        rep.write(f"Mapping: {'wfmash --approx-map' if mapped else 'none (grouping by sequence names)'}\n\n")
        for r in rows:
            rep.write(f"{r[0]}: {r[2]} sequences, {r[3]:,} bp, {r[1]} genomes\n")
        rep.write("\n")
        for r in rows:
            members = by_group[r[0]]
            count = collections.Counter(genome_of[n] for n in members)
            for gen in genomes:
                if r[0] != "unplaced" and count[gen] == 0:
                    rep.write(f"WARNING: no sequence of {gen} in group {r[0]}\n")
                elif r[0] != "unplaced" and count[gen] > 1:
                    rep.write(f"NOTE: {count[gen]} sequences of {gen} in group {r[0]}\n")
        low = [n for n in query_seqs if assignment[n][1] in ("mapping_low", "none", "name")]
        for n in low:
            g, method, best, frac = assignment[n]
            rep.write(f"NOTE: {n} ({length[n]:,} bp) -> {g} [{method}, best={best}, mapped fraction={frac}]\n")
    sys.stderr.write(open(os.path.join(args.outdir, "partition_report.txt")).read())


if __name__ == "__main__":
    main()
