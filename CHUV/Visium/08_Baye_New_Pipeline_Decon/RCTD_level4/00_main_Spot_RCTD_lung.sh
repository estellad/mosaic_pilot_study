#!/bin/bash
for i in `seq 1 5`;
do 
    file="/work/PRTNR/CHUV/DIR/rgottar1/spatial/env/ydong/mosaic_pilot_study/out_08_Spot_Decon_RCTD_lung_${i}"
    sbatch --output=$file /work/PRTNR/CHUV/DIR/rgottar1/spatial/env/ydong/mosaic_pilot_study/CHUV/Visium/08_Baye_New_Pipeline_Decon/RCTD_level4/00_Spot_Decon_Level4_RCTD_lung.sh  $i
done