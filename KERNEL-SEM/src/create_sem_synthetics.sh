#!/bin/bash
#
# Test 1 - Plane Wave Example
#

echo "running Test 1: `date`"
currentdir=`pwd`

NPROC=1

DIP_ANGLE_DEGREES=0.0

LAB_WIDTH=0.0

DELTA_T=0.02
T_END=300  #300 usually good

#must make input files
python << END
import sys
sys.path.append('src/')
import write_input_files as w
import create_tomo_file as tomo
H=1500000.0
W=3000000.0
N_ELEM_VERT=751
N_ELEM_HORIZ=1501
params={
"N_STATIONS" : 150,
"HEIGHT" : H,
"WIDTH" : W,
"DELTA_T" : $DELTA_T,
"T_END" : $T_END,
"DEPTH_TO_LAB" : 1000000000.0,
"DEPTH_TO_MOHO" :     00000.0,
"VERTICAL_SMOOTHING_LAMBDA" : 0000.0,
"HORIZONTAL_SMOOTHING_LAMBDA" :0000.0,
"SCAT_DEPTH" : $1,
"SCAT_STRENGTH" : $2,
"SCAT_RADIUS" : $3,
"ILLUMINATION_DIRECTION" : 1,
"GAP_WIDTH" : 0.0,
"ANGLE_SOURCE" : $4
}
params["NTSTEP"]=round( params["T_END"] / params["DELTA_T"] )
tomo.main(**params)
w.write_Par_file($NPROC,H,N_ELEM_VERT,N_ELEM_HORIZ,$LAB_WIDTH,W,**params)
w.write_SOURCE(H,**params)
w.write_interfaces(H,N_ELEM_VERT,W,$DIP_ANGLE_DEGREES,$LAB_WIDTH)
END

echo
echo "   setting up example..."
echo

mkdir -p OUTPUT_FILES
mkdir -p DATA

# sets up local DATA/ directory
cd DATA/
rm -f Par_file SOURCE
ln -s ../Par_file.in Par_file
ln -s ../SOURCE.in SOURCE
cd ../

# cleans output files
rm -rf OUTPUT_FILES/*

cd $currentdir

# links executables
rm -f xmeshfem2D xspecfem2D
XPATH='/users/nmancine/PROG/SPECFEM2D/specfem2d/bin'
ln -s $XPATH/xmeshfem2D
ln -s $XPATH/xspecfem2D

# stores setup
cp DATA/Par_file OUTPUT_FILES/
cp DATA/SOURCE OUTPUT_FILES/
cp interfaces.dat OUTPUT_FILES/

# Get the number of processors
NPROC=`grep ^NPROC DATA/Par_file | cut -d = -f 2 | cut -d \# -f 1 | tr -d ' '`

# runs database generation
echo
echo "  running mesher..."
echo
./xmeshfem2D

if [ "$NPROC" -eq 1 ]; then # This is a serial simulation
  echo
  echo " Running solver..."
  echo
  ./xspecfem2D
else # This is a MPI simulation
  echo
  echo " Running solver on $NPROC processors..."
  echo
  mpirun -np $NPROC ./xspecfem2D
fi

# stores output
cp DATA/*SOURCE* DATA/*STATIONS* OUTPUT_FILES

echo
echo "see results in directory: OUTPUT_FILES/"
echo
echo "done"

date
