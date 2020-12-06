#!/bin/bash

for inst in ../instances/*.txt; do
  qsub -l h_rt=00:16:00 -r Y -e out.txt $1 $(basename "$inst")
done;