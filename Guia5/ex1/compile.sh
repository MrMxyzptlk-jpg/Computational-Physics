gfortran ./modulos/precision.f90\
        ./modulos/mzranmod.f90\
        ./modulos/formats.f90\
        ./modulos/constantes.f90\
        ./modulos/funciones.f90\
        ./modulos/subrutinas.f90\
        ./modulos/linkedLists.f90\
        ./modulos/parsing.f90\
        ./modulos/writing2files.f90\
        ex1.f90 -o ex1.exe -ffpe-trap=invalid,overflow,zero \
         -O3 -march=native -ftree-vectorize -ftree-vectorizer-verbose=2 \
         -O -Wall -fcheck=all -g -fbacktrace \
         -fopenmp -O2

echo
