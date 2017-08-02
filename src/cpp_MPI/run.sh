#!/bin/bash

rm *.txt

make clean

make -j 7

echo "compilation & linkin done"

echo "Running model.out"

#export OMP_NUM_THREADS=2

mpirun -n 4 ./model.out > output_test.txt

echo "DONE"
