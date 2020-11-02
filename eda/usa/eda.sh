#!/usr/bin/env bash
CSV=../../data/usa_ecoli_shigella_2000_on.tsv
OUTDIR=../../results/
./map.py --csv $CSV --outdir $OUTDIR
