#!/bin/bash
mpicc -o $1.run $1.c
mpirun --oversubscribe -n 8 ./$1.run
