gfortran ./modulos/precision.f90\
        ./modulos/formats.f90\
        ./modulos/mzranmod_threadsafe.f90\
        ./modulos/subrutinas.f90\
        ./modulos/parsing.f90\
        ex1.f90 -o ex1.exe  \
        -O3 -march=native -ftree-vectorize -ftree-vectorizer-verbose=2 \
        -ffpe-trap=invalid,overflow,zero -O -Wall -fcheck=all -g -fbacktrace\
        -fopenmp -O2

#./ex1.exe

echo
