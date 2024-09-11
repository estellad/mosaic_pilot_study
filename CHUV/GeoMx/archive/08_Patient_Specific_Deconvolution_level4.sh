#!/bin/bash
#SBATCH --output out_Geo_Pt_Spec_Level4_Decon
#SBATCH --partition cpu
#SBATCH --time=72:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=100gb
module load singularityce
singularity exec --bind /scratch,/users,/work /work/PRTNR/CHUV/DIR/rgottar1/spatial/containers/estelladong/rserver_sk_master.sif Rscript --vanilla /work/PRTNR/CHUV/DIR/rgottar1/spatial/env/ydong/mosaic_pilot_study/CHUV/GeoMx/08_Patient_Specific_Deconvolution_level4.R