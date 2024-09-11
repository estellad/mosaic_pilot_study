#!/bin/bash
for i in `seq 1 6`;
do 
    file="/work/PRTNR/CHUV/DIR/rgottar1/spatial/env/ydong/mosaic_pilot_study/out_08_Spot_Decon_RCTD_dlbcl_${i}"
    sbatch --output=$file /work/PRTNR/CHUV/DIR/rgottar1/spatial/env/ydong/mosaic_pilot_study/CHUV/Visium/08_Baye_New_Pipeline_Decon/RCTD/00_Spot_Decon_Level1_5_RCTD_dlbcl.sh $i
done