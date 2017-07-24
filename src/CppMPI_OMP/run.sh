#!/bin/bash

export OMP_NUM_THREADS=2

mpiexec.hydra -n 2 ./model.out > output.txt
