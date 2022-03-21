#!/bin/bash
echo "Problem Size;Sequential;Parallel" > ../../plotting/$1.csv
for i in {1..15}
do
    ./$1.run $i >> ../../plotting/$1.csv
done