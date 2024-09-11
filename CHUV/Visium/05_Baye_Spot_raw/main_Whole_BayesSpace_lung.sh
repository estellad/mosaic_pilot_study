#!/bin/bash
for i in `seq 1 5`;
do 
    file="/work/PRTNR/CHUV/DIR/rgottar1/spatial/env/ydong/Owkin_Pilot/Code/Visium/Manuscript/05_Baye_Spot_raw/out_05_Baye_Spot_raw_Whole_BayesSpace_lung_${i}"
    sbatch --output=$file /work/PRTNR/CHUV/DIR/rgottar1/spatial/env/ydong/Owkin_Pilot/Code/Visium/Manuscript/05_Baye_Spot_raw/Whole_BayesSpace_lung.sh $i
done