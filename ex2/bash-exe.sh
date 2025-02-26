gfortran ./modulos/ISO.f90\
        ./modulos/precision.f90\
        ./modulos/funciones.f90\
        ./modulos/subrutinas.f90\
        ex2.f90 -o ex2.exe 
./ex2.exe

gnuplot plot.gp

echo