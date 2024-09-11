#!/bin/bash
#SBATCH --output out_level4_RCTD_B3_2
#SBATCH --partition cpu
#SBATCH --time=8:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=20gb
module load singularityce
singularity exec --bind /scratch,/users,/work /work/PRTNR/CHUV/DIR/rgottar1/spatial/containers/estelladong/rserver_sk_master_noglmGamPoi_mosaicpaper.sif Rscript --vanilla /work/PRTNR/CHUV/DIR/rgottar1/spatial/env/ydong/mosaic_pilot_study/CHUV/Visium/08_Baye_New_Pipeline_Decon/RCTD/00_Spot_Decon_Level4_all_RCTD.R