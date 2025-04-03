gfortran ./modulos/precision.f90\
        ./modulos/formats.f90\
        ./modulos/constantes.f90\
        ./modulos/funciones.f90\
        ./modulos/subrutinas.f90\
        ex4.f90 -o ex4.exe -ffpe-trap=invalid,overflow,zero -lfftw3 -lfftw3_threads -O -Wall -fcheck=all -g -fbacktrace
./ex4.exe

#gnuplot plot.gp

echo
