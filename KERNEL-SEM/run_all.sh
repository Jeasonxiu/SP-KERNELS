#!/bin/bash
#
for SCAT_STRENGTH in 0.01; do
for SCAT_RADIUS in 2500 5000 7500; do
for ANGLE_SOURCE in 15 20 25; do
for SCAT_DEPTH in `seq 000 5000 300000`; do
DIRNAME=$SCAT_DEPTH-$SCAT_STRENGTH-$SCAT_RADIUS-$ANGLE_SOURCE-DBLPERIOD
sbatch run.sh $DIRNAME $SCAT_DEPTH $SCAT_STRENGTH $SCAT_RADIUS $ANGLE_SOURCE
exit 0
done 
done
done
done
#