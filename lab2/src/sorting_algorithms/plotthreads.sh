#!/bin/bash
echo "Problem Size;2;4;8;16" > ../../plotting/$1.csv
for i in {1..20}
do
    echo -n $i >> ../../plotting/$1.csv
    make -B OMP_NUM_THREADS=2 > /dev/null
    ./$1.run $i >> ../../plotting/$1.csv
    make -B OMP_NUM_THREADS=4 > /dev/null
    ./$1.run $i >> ../../plotting/$1.csv
    make -B OMP_NUM_THREADS=8 > /dev/null
    ./$1.run $i >> ../../plotting/$1.csv
    make -B OMP_NUM_THREADS=16 > /dev/null
    ./$1.run $i >> ../../plotting/$1.csv
    echo ""  >> ../../plotting/$1.csv 
done