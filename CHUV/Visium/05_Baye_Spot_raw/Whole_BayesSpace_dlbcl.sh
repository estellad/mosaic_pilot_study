#!/bin/bash
#SBATCH --partition cpu
#SBATCH --time=72:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=100gb
module load singularity
singularity exec --bind /scratch,/users,/work /work/PRTNR/CHUV/DIR/rgottar1/spatial/containers/estelladong/rserver_sk_master.sif Rscript --vanilla /work/PRTNR/CHUV/DIR/rgottar1/spatial/env/ydong/Owkin_Pilot/Code/Visium/Manuscript/05_Baye_Spot_raw/Whole_BayesSpace_dlbcl.R $1