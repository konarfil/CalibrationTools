#!/bin/sh
#SBATCH --mem      8192M
#SBATCH --licenses sps
#SBATCH --time     12:00:00
#SBATCH --ntasks   1

./exctract_OM_spectra -i $RED_PATH/snemo_run-%RUN_NUMBER_red-v2.data.gz -o %OUTPUT_FOLDER -r %RUN_NUMBER -n %N_EVENTS
root '../fit_spectra.cpp("%OUTPUT_FOLDER", "%FIT_FOLDER")'
