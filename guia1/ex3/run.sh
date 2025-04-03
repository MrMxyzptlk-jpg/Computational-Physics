gfortran ./modulos/precision.f90\
        ./modulos/formats.f90\
        ./modulos/constantes.f90\
        ./modulos/funciones.f90\
        ./modulos/subrutinas.f90\
        ex3.f90 -o ex3.exe -O3 -march=native -ftree-vectorize -ftree-vectorizer-verbose=2 \
         -lfftw3 -lfftw3_threads -ffpe-trap=invalid,overflow,zero -O -Wall -fcheck=all -g -fbacktrace  
./ex3.exe

echo
