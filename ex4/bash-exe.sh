gfortran ./modulos/ISO.f90\
        ./modulos/precision.f90\
        ./modulos/funciones.f90\
        ./modulos/subrutinas.f90\
        ex4.f90 -o ex4.exe 
./ex4.exe
gnuplot plot.gp

echo