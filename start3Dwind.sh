#!/bin/bash



export outputfolder=${2}

fpath='./'$outputfolder'/'temp_packet*.dat
prevfile0='./'$outputfolder'/'packet000.dat
prevfile='./'$outputfolder'/'packet*.dat

if [ -e $prevfile0 ]; then
 rm $prevfile
fi

mpirun.mpich -np $1 ./main.run

rm $fpath

cp input.dat $outputfolder

inputmodelfile=`cat input.dat | grep inputmodelFile | cut -d= -f 2`
cp $inputmodelFile $outputfolder/

python3.6 spec.py $outputfolder
