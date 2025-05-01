gfortran ./modulos/precision.f90\
        ./modulos/formats.f90\
        ./modulos/subrutinas.f90\
        ./modulos/funciones.f90\
        ./modulos/mzranmod.f90\
        ex3.f90 -o ex3.exe -ffpe-trap=invalid,overflow,zero
#./ex3.exe

echo
