#!/bin/bash
mpicc -o $1.run $1.c
mpirun -n 4 ./$1.run
