#!/usr/bin/env bash
CSV=../../data/ECDC_surveillance_data_Antimicrobial_resistance.csv
OUTDIR=../../results/
./corr.py --csv $CSV --outdir $OUTDIR
