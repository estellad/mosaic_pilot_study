#!/bin/bash
#SBATCH --partition cpu
#SBATCH --time=10:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=50gb
module load singularityce
singularity exec --bind /scratch,/users,/work /work/PRTNR/CHUV/DIR/rgottar1/spatial/containers/estelladong/rserver_sk_master_noglmGamPoi_mosaicpaper.sif Rscript --vanilla /work/PRTNR/CHUV/DIR/rgottar1/spatial/env/ydong/mosaic_pilot_study/CHUV/Visium/08_Baye_New_Pipeline_Decon/RCTD/00_Spot_Decon_Level1_5_RCTD_lung.R $1