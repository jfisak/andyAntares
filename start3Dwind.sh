#!/bin/bash



export OUTPUTFO=${2}

FOLDERPATH='./'${2}'/'temp_packet*.dat
PREVFILE0='./'${2}'/'packet000.dat
PREVFILE='./'${2}'/'packet*.dat

if [ -e $PREVFILE0 ]; then
 rm $PREVFILE
fi

mpirun.mpich -np $1 ./main.run

rm $FOLDERPATH

