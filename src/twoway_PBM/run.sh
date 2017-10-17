#!/bin/bash

rm -f *.txt

make clean

make -j 7

echo "compilation & linkin done"

echo "Running model.out"

#export OMP_NUM_THREADS=2

mpirun -n 1 ./model.out 128 200 > output_test.txt

echo "DONE"

#matlab -nodisplay -nosplash -r "d50andRatioplt;quit";
#matlab -nodisplay -nosplash  exit;

#echo "Plots saved in csvDump folder"
