#!/usr/bin/env bash
# panexplorer_auto.sh - "auto" mode of PanExplorer2
#
# 1. runs Snakemake_files/Snakefile_route: prepares the genomes, computes the ANI between all genomes
#    (FastANI for prokaryotes, Mash for eukaryotes) and chooses the pangenome method;
# 2. runs the pangenome workflow that was chosen (PGGB if the minimum ANI >= threshold, otherwise
#    PanACoTA / Roary for prokaryotes or OrthoFinder for eukaryotes).
# The genomes already prepared in outputs/genomes are reused by the second workflow.
#
# Usage (from the directory containing config.yaml):
#   bash $PANEX_PATH/panexplorer_auto.sh [--input-type prokaryote|eukaryote] [--threshold 95]
#        [--force auto|graph|clustering] [--clustering-tool panacota|roary|orthofinder]
#        [--cores N] [--config key=value ...]
# The decision is written in outputs/route/decision.txt (and decision.tsv, ani_tree.nwk, ani_pairs.tsv).

set -euo pipefail

PANEX_PATH="${PANEX_PATH:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
SNAKE_DIR="$PANEX_PATH/Snakemake_files"

# empty = value of config.yaml (pangenome_method, ani_threshold, clustering_tool), else the defaults
INPUT_TYPE="prokaryote"
THRESHOLD=""
FORCE=""
CLUSTERING_TOOL=""
CORES="1"
EXTRA_CONFIG=()

usage() { sed -n '2,15p' "$0" | sed 's/^# \{0,1\}//'; }

while [ $# -gt 0 ]; do
    case "$1" in
        --input-type)      INPUT_TYPE="$2"; shift 2 ;;
        --threshold)       THRESHOLD="$2"; shift 2 ;;
        --force)           FORCE="$2"; shift 2 ;;
        --clustering-tool) CLUSTERING_TOOL="$2"; shift 2 ;;
        --cores)           CORES="$2"; shift 2 ;;
        --config)          shift
                           while [ $# -gt 0 ] && [ "${1#--}" = "$1" ]; do EXTRA_CONFIG+=("$1"); shift; done ;;
        -h|--help)         usage; exit 0 ;;
        *) echo "Unknown option: $1" >&2; usage >&2; exit 1 ;;
    esac
done

case "$INPUT_TYPE" in prokaryote|eukaryote) ;; *) echo "--input-type must be prokaryote or eukaryote" >&2; exit 1 ;; esac
case "$FORCE" in ""|auto|graph|clustering) ;; *) echo "--force must be auto, graph or clustering" >&2; exit 1 ;; esac
ROUTE_CONFIG=(input_type="$INPUT_TYPE")
[ -n "$THRESHOLD" ] && ROUTE_CONFIG+=(ani_threshold="$THRESHOLD")
[ -n "$FORCE" ] && ROUTE_CONFIG+=(pangenome_method="$FORCE")
[ -n "$CLUSTERING_TOOL" ] && ROUTE_CONFIG+=(clustering_tool="$CLUSTERING_TOOL")
if [ ! -f config.yaml ]; then echo "config.yaml not found in $(pwd)" >&2; exit 1; fi

# do not rerun the preparation steps only because the Snakefile or its parameters differ
RERUN_OPT=()
if snakemake --help 2>/dev/null | grep -q -- "--rerun-triggers"; then RERUN_OPT=(--rerun-triggers mtime); fi

echo "=== PanExplorer2 auto mode: ANI-based choice of the pangenome method ==="
echo "options: ${ROUTE_CONFIG[*]}"

snakemake --cores "$CORES" -s "$SNAKE_DIR/Snakefile_route" \
    --config "${ROUTE_CONFIG[@]}" ${EXTRA_CONFIG[@]+"${EXTRA_CONFIG[@]}"}

METHOD=$(tr -d '[:space:]' < outputs/route/method.txt)
cat outputs/route/decision.txt

case "$INPUT_TYPE:$METHOD" in
    prokaryote:pggb)        SNAKEFILE="Snakefile_wget_pggb_heatmap_upset_COG" ;;
    prokaryote:panacota)    SNAKEFILE="Snakefile_wget_panacota_heatmap_upset_COG" ;;
    prokaryote:roary)       SNAKEFILE="Snakefile_wget_roary_heatmap_upset_COG" ;;
    eukaryote:pggb)         SNAKEFILE="Snakefile_pggb_heatmap_upset_COG" ;;
    eukaryote:orthofinder)  SNAKEFILE="Snakefile_orthofinder_heatmap_upset" ;;
    *) echo "No workflow for input type '$INPUT_TYPE' and method '$METHOD'" >&2; exit 1 ;;
esac

echo "=== Running $SNAKEFILE ==="
snakemake --cores "$CORES" ${RERUN_OPT[@]+"${RERUN_OPT[@]}"} -s "$SNAKE_DIR/$SNAKEFILE" \
    ${EXTRA_CONFIG[@]+"--config"} ${EXTRA_CONFIG[@]+"${EXTRA_CONFIG[@]}"}
