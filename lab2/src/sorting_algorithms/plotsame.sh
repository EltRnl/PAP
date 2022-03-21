#!/bin/bash
echo "Problem Size;Sequential;Parallel" > ../../plotting/$1.csv
for i in {1..20}
do
    ./$1.run $i >> ../../plotting/$1.csv
done
