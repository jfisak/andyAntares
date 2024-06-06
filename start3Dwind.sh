#!/bin/bash

# basic variables definition
program='mpirun.mpich'
pythonVersion='python3.11'

outputfolder=${2}

# the name will be created based on date
echo $outputfolder
if [[ $outputfolder = '' ]]; then
 datum=`date`
 year=`echo $datum | cut -d' ' -f6`
 month=`echo $datum | cut -d' ' -f3`
 nday=`echo $datum | cut -d' ' -f2`
 hour=`echo $datum | cut -d' ' -f4 | cut -d ':' -f1`
 minute=`echo $datum | cut -d' ' -f4 | cut -d ':' -f2`
 second=`echo $datum | cut -d' ' -f4 | cut -d ':' -f3`
 # create a name of an output folder
 outputfolder='3dwind'$year$month$nday$hour$minute$second
fi
export outputfolder

# test of exsistence of the folder
if [ ! -e $outputfolder ]; then
 mkdir $outputfolder
fi

FOLDERPATH='./'$outputfolder'/'temp_packet*.dat
PREVFILE0='./'$outputfolder'/'packet000.dat
PREVFILE='./'$outputfolder'/'packet*.dat

if [ -e $PREVFILE0 ]; then
 rm $PREVFILE
fi

$program -np $1 ./main.run

rm $FOLDERPATH

cp input.dat $outputfolder

$pythonVersion spec2.py $outputfolder
