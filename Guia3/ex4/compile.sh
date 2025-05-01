gfortran ./modulos/precision.f90\
        ./modulos/formats.f90\
        ./modulos/constantes.f90\
        ./modulos/mzranmod.f90\
        ./modulos/mzranmod_threadsafe.f90\
        ./modulos/mtmod.f90\
        ./modulos/funciones.f90\
        ./modulos/subrutinas.f90\
        ex4.f90 -o ex4.exe  -O3 -march=native -ftree-vectorize -ftree-vectorizer-verbose=2 \
                -ffpe-trap=invalid,overflow,zero -O -Wall -fcheck=all -g -fbacktrace \
                -fopenmp -O2
#./ex4.exe

echo
