#!/bin/bash
cd ..
python3.7 solve.py --alg vnd --inst_file $1 --mh_ttime $((15*60)) --out_file res/vnd.csv