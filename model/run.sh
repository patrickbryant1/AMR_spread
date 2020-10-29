#!/usr/bin/env bash
CSV=../data/ECDC_surveillance_data_Antimicrobial_resistance.csv
OUTDIR=../results/
./calc_prob.py --csv $CSV --outdir $OUTDIR
