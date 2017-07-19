#!/bin/bash

rm output.txt

make clean

make -j 6

echo "compilation & linkin done"

echo "Running model.out"

export OMP_NUM_THREADS=2

mpirun -n 4 ./model_256_555.out > output_test.txt

echo "DONE"
