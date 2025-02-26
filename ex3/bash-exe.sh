gfortran ./modulos/ISO.f90\
        ./modulos/precision.f90\
        ./modulos/funciones.f90\
        ./modulos/subrutinas.f90\
        ex3.f90 -o ex3.exe 
./ex3.exe

gnuplot plot.gp

echo