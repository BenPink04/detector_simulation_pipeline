#!/bin/bash

for i in {1..10}
do
    SEED=$((1000000 + $i))
    OUTFILE="out_seed_Scen4_${i}.root"
    # Run Geant4 with macro file and seed
    ./jlabk run.mac $SEED
    # Move or rename the output file
    mv out.root $OUTFILE
done