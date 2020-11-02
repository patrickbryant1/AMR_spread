#!/usr/bin/env bash
TSV=../../data/isolates.tsv #usa_ecoli_shigella_2000_on.tsv
OUTDIR=../../results/usa/
./map.py --tsv $TSV --outdir $OUTDIR
