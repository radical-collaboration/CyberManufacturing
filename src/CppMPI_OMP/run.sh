#!/bin/bash

export OMP_NUM_THREADS=2

mpiexec.hydra -n 3 ./model.out > output.txt
