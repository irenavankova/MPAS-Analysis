#!/bin/bash
#SBATCH  --job-name=mpas_analysis
#SBATCH  --nodes=2
#SBATCH  --output=mpas_analysis.o%j
#SBATCH  --exclusive
#SBATCH  --time=1:00:00
#SBATCH  --qos=standard
#SBATCH  --partition=standard


source /users/vankova/miniconda3/etc/profile.d/conda.sh
conda activate mpas_dev

export HDF5_USE_FILE_LOCKING=FALSE

mpas_analysis --verbose -m chicoma-cpu /users/vankova/mpas-analysis/fris_map2/cfg_files/sorrmGcases/chic_e3sm_test_20_29.cfg





