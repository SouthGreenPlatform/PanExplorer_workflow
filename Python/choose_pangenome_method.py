#!/usr/bin/env python3
"""
choose_pangenome_method.py - decide between the graph-based (PGGB) and the clustering-based pangenome
workflows from the genome-to-genome similarity (ANI), and write a small report + ANI-based NJ tree.

Rationale: the gene presence/absence matrix derived from a PGGB graph reproduces gene clustering
(PanACoTA/OrthoFinder) within a species, but no longer links orthologs between species: below ~95% ANI
most homologous regions are not aligned into shared segments. Within a species (min ANI >= threshold)
the graph mode is used, otherwise the clustering mode.

Input (one of):
  --fastani FILE   fastANI output (query, reference, ANI, fragments, total). Pairs that FastANI does not
                   report (ANI below ~80%) are considered as divergent (ANI < 80).
  --mash FILE      "mash triangle" output (distances); ANI is estimated as 100 * (1 - distance).
Options:
  --strains FILE   two-column file (id <TAB> displayed name), e.g. outputs/genomes/strains.txt
  --threshold 95   minimum ANI (%) between all genomes to use the graph mode
  --force auto|graph|clustering
  --graph-tool pggb   --clustering-tool panacota|roary|orthofinder
  --outdir outputs/route
Outputs (in --outdir):
  method.txt        one word: the tool to run (e.g. pggb, panacota, orthofinder)
  decision.tsv      key <TAB> value (mode, tool, min/mean ANI, threshold, closest/most distant pairs...)
  decision.txt      human-readable explanation (also printed)
  ani_pairs.tsv     genome1 <TAB> genome2 <TAB> ANI (symmetrised; NA if not reported)
  ani_tree.nwk      neighbour-joining tree on (100 - ANI) distances
"""
import argparse
import math
import os
import re
import sys

MISSING_ANI = 80.0  # value used for pairs not reported by FastANI (they are below its detection limit)


def genome_id(path):
    """outputs/genomes/GCF_000019585.2.fasta -> GCF_000019585.2"""
    base = os.path.basename(path.strip())
    return re.sub(r"\.(fasta|fna|fa|fas)(\.gz)?$", "", base)


def read_strains(path):
    names = {}
    if path and os.path.exists(path):
        with open(path) as fh:
            for line in fh:
                parts = line.rstrip("\n").split("\t")
                if len(parts) >= 2 and parts[0]:
                    names[parts[0]] = parts[1]
    return names


def read_fastani(path):
    genomes, ani = set(), {}
    with open(path) as fh:
        for line in fh:
            parts = line.split()
            if len(parts) < 3:
                continue
            a, b, v = genome_id(parts[0]), genome_id(parts[1]), float(parts[2])
            genomes.update((a, b))
            if a != b:
                ani.setdefault(frozenset((a, b)), []).append(v)
    return sorted(genomes), {k: sum(v) / len(v) for k, v in ani.items()}


def read_mash_triangle(path):
    with open(path) as fh:
        lines = [l.rstrip("\n") for l in fh if l.strip()]
    genomes, ani = [], {}
    for line in lines[1:]:  # first line = number of sequences
        parts = line.split("\t")
        g = genome_id(parts[0])
        for j, d in enumerate(parts[1:]):
            ani[frozenset((g, genomes[j]))] = 100.0 * (1.0 - float(d))
        genomes.append(g)
    return genomes, ani


def nj_newick(names, dist):
    """Neighbour-joining (Saitou & Nei 1987) on a full distance matrix (list of lists)."""
    nodes = [n.replace(":", "_").replace(",", "_").replace("(", "_").replace(")", "_").replace(" ", "_")
             for n in names]
    D = [row[:] for row in dist]
    if len(nodes) == 1:
        return nodes[0] + ";"
    while len(nodes) > 2:
        n = len(nodes)
        r = [sum(D[i]) for i in range(n)]
        best, bi, bj = None, 0, 1
        for i in range(n):
            for j in range(i + 1, n):
                q = (n - 2) * D[i][j] - r[i] - r[j]
                if best is None or q < best:
                    best, bi, bj = q, i, j
        li = max(0.0, 0.5 * D[bi][bj] + (r[bi] - r[bj]) / (2 * (n - 2)))
        lj = max(0.0, D[bi][bj] - li)
        new = f"({nodes[bi]}:{li:.5f},{nodes[bj]}:{lj:.5f})"
        newrow = [0.5 * (D[bi][k] + D[bj][k] - D[bi][bj]) for k in range(n)]
        keep = [k for k in range(n) if k not in (bi, bj)]
        D = [[D[a][b] for b in keep] + [newrow[a]] for a in keep] + [[newrow[b] for b in keep] + [0.0]]
        nodes = [nodes[k] for k in keep] + [new]
    return f"({nodes[0]}:{D[0][1] / 2:.5f},{nodes[1]}:{D[0][1] / 2:.5f});"


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    src = ap.add_mutually_exclusive_group(required=True)
    src.add_argument("--fastani")
    src.add_argument("--mash")
    ap.add_argument("--strains")
    ap.add_argument("--threshold", type=float, default=95.0)
    ap.add_argument("--force", choices=["auto", "graph", "clustering"], default="auto")
    ap.add_argument("--graph-tool", default="pggb")
    ap.add_argument("--clustering-tool", default="panacota")
    ap.add_argument("--outdir", default="outputs/route")
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    if args.fastani:
        genomes, ani = read_fastani(args.fastani)
        method_ani = "FastANI"
    else:
        genomes, ani = read_mash_triangle(args.mash)
        method_ani = "Mash (ANI = 100 x (1 - distance))"
    names = read_strains(args.strains)
    label = lambda g: names.get(g, g)

    if len(genomes) < 2:
        sys.exit("[choose_pangenome_method] fewer than 2 genomes found in the ANI file")

    pairs, missing = [], 0
    for i, a in enumerate(genomes):
        for b in genomes[i + 1:]:
            v = ani.get(frozenset((a, b)))
            if v is None:
                missing += 1
            pairs.append((a, b, v))
    values = [v if v is not None else MISSING_ANI for _, _, v in pairs]
    min_ani, mean_ani, max_ani = min(values), sum(values) / len(values), max(values)
    closest = max(pairs, key=lambda p: p[2] if p[2] is not None else MISSING_ANI)
    farthest = min(pairs, key=lambda p: p[2] if p[2] is not None else MISSING_ANI)

    within_species = min_ani >= args.threshold
    if args.force == "graph":
        mode, reason = "graph", "forced by the user"
    elif args.force == "clustering":
        mode, reason = "clustering", "forced by the user"
    else:
        mode = "graph" if within_species else "clustering"
        reason = (f"minimum ANI {min_ani:.2f}% >= {args.threshold:g}%: genomes of the same species, "
                  f"graph-based gene PAV is reliable" if within_species else
                  f"minimum ANI {min_ani:.2f}% < {args.threshold:g}%: genomes from different species, "
                  f"a graph would not link orthologs between species -> gene clustering")
    tool = args.graph_tool if mode == "graph" else args.clustering_tool
    warning = ""
    if mode == "graph" and not within_species:
        warning = (f"WARNING: graph mode forced although the minimum ANI ({min_ani:.2f}%) is below "
                   f"{args.threshold:g}%; orthologs of distant genomes will mostly appear as distinct clusters.")

    with open(os.path.join(args.outdir, "method.txt"), "w") as fh:
        fh.write(tool + "\n")
    with open(os.path.join(args.outdir, "ani_pairs.tsv"), "w") as fh:
        fh.write("genome1\tgenome2\tANI\n")
        for a, b, v in pairs:
            fh.write(f"{label(a)}\t{label(b)}\t{'NA' if v is None else f'{v:.4f}'}\n")
    decision = {
        "mode": mode, "tool": tool, "reason": reason, "ani_method": method_ani,
        "threshold": f"{args.threshold:g}", "n_genomes": len(genomes), "n_pairs": len(pairs),
        "pairs_not_reported": missing, "min_ANI": f"{min_ani:.4f}", "mean_ANI": f"{mean_ani:.4f}",
        "max_ANI": f"{max_ani:.4f}",
        "most_distant_pair": f"{label(farthest[0])} / {label(farthest[1])}",
        "closest_pair": f"{label(closest[0])} / {label(closest[1])}",
        "forced": args.force, "warning": warning,
    }
    with open(os.path.join(args.outdir, "decision.tsv"), "w") as fh:
        for k, v in decision.items():
            fh.write(f"{k}\t{v}\n")

    # NJ tree on 100 - ANI distances (missing pairs set to 100 - MISSING_ANI)
    idx = {g: i for i, g in enumerate(genomes)}
    n = len(genomes)
    D = [[0.0] * n for _ in range(n)]
    for a, b, v in pairs:
        d = 100.0 - (v if v is not None else MISSING_ANI)
        D[idx[a]][idx[b]] = D[idx[b]][idx[a]] = d
    with open(os.path.join(args.outdir, "ani_tree.nwk"), "w") as fh:
        fh.write(nj_newick([label(g) for g in genomes], D) + "\n")

    text = [
        "PanExplorer2 - choice of the pangenome method",
        f"  genomes: {len(genomes)}   pairs: {len(pairs)}   ANI computed with {method_ani}",
        f"  ANI min / mean / max: {min_ani:.2f}% / {mean_ani:.2f}% / {max_ani:.2f}%"
        + (f"   ({missing} pairs below the FastANI detection limit, counted as {MISSING_ANI:g}%)" if missing else ""),
        f"  most distant pair: {decision['most_distant_pair']}",
        f"  threshold: {args.threshold:g}%   requested: {args.force}",
        f"  -> mode: {mode}   tool: {tool}",
        f"  reason: {reason}",
    ]
    if warning:
        text.append("  " + warning)
    with open(os.path.join(args.outdir, "decision.txt"), "w") as fh:
        fh.write("\n".join(text) + "\n")
    print("\n".join(text))


if __name__ == "__main__":
    main()
