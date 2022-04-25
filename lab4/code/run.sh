#!/bin/bash

case $1 in
    lbm)
        make
        mpirun -np 4 ./lbm -e $2
        ;;
    com)
        shift
        mpirun -np 4 ./check_comm -e $@
        ;;
    gif)
        ./gen_animate_gif.sh output.raw output.gif
        ;;
    gifm)
        ./gen_animate_gif.sh $2.raw $3.gif
        ;;
    verif)
        ./display --gnuplot $2.raw 3 | md5sum -c ref_ex0.md5
        ;;
    clean)
        make clean
        rm *.gif *.raw
        ;;
    *)
        echo "Usage :"
        printf "\t lbm <nb_exo>\n\t ↳ makes and run lbm with the exercice number given as argument\n"
        printf "\t com [args]\n\t ↳ runs the communcation debuging with all the options given in argument (if any)\n"
        printf "\t gif\n\t ↳ transform 'output.raw' into 'output.gif'\n"
        printf "\t gifm <input>.raw <output>.gif\n\t ↳ same as previous, but will take the names 'input' and 'output' given to load the raw and save the gif\n\t\t note : you don't need to add the extensions, it's already doing it\n"
        printf "\t verif <input>.raw\n\t ↳ checks if the output is the same as exercice 1\n"

        ;;
esac