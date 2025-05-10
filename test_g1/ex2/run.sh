# This might not be the proper compilation procedure but it works. See the other "run.sh" files in the other exercises of "/guia1" for better examples.
gfortran ./modulos/precision.f90\
        ./modulos/funciones.f90\
        ./modulos/subrutinas.f90\
        ex2.f90 -o ex2.exe 

gfortran -c  ./modulos/funciones.f90\
        ./modulos/subrutinas.f90\

gfortran ex2.f90 funciones.o subrutinas.o -o ex2.exe -lfftw3 -lfftw3_threads -ffpe-trap=invalid,overflow,zero


./ex2.exe             # Executable file
