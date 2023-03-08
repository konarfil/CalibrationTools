#!/bin/sh

echo "Choose run number:"
read RUN_NUMBER
echo "					   "

echo "Choose how many events to read (-1 to read all events of the run):"
read N_EVENTS
echo "					   "

mkdir -p ../runs/run_$RUN_NUMBER/OM_histos
mkdir -p ../runs/run_$RUN_NUMBER/fits

PARENT_DIR=$( cd ..;pwd)

sed -e "s|%RUN_NUMBER|$RUN_NUMBER|g" \
    -e "s|%N_EVENTS|$N_EVENTS|g" \
    -e "s|%OUTPUT_FOLDER|$PARENT_DIR/runs/run_$RUN_NUMBER/OM_histos/|g" \
    -e "s|%FIT_FOLDER|$PARENT_DIR/runs/run_$RUN_NUMBER/fits/|g" \
    ./temp_sub.sh > $PARENT_DIR/runs/run_$RUN_NUMBER/sub.sh
    
sbatch -o $PARENT_DIR/runs/run_$RUN_NUMBER/OUT.log -e $PARENT_DIR/runs/run_$RUN_NUMBER/ERR.log $PARENT_DIR/runs/run_$RUN_NUMBER/sub.sh
