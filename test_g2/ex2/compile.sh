gfortran ./modulos/precision.f90\
        ./modulos/formats.f90\
        ./modulos/constantes.f90\
        ./modulos/funciones.f90\
        ./modulos/subrutinas.f90\
        ex2.f90 -o ex2.exe -ffpe-trap=invalid,overflow,zero -O -Wall -fcheck=all -g -fbacktrace

# ./ex2.exe # Executable file

echo
