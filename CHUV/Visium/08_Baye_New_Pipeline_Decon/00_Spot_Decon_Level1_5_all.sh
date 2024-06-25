#!/bin/bash
#SBATCH --output out_08_Spot_DLBCL_Decon_new_pipelines
#SBATCH --partition cpu
#SBATCH --time=72:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=100gb
module load singularity
singularity exec --bind /scratch,/users,/work /work/PRTNR/CHUV/DIR/rgottar1/spatial/containers/estelladong/rserver_sk_master.sif Rscript --vanilla /work/PRTNR/CHUV/DIR/rgottar1/spatial/env/ydong/Owkin_Pilot/Code/Visium/Manuscript/08_Baye_New_Pipeline_Decon/00_Spot_Decon_Level1_5_all.R