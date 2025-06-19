gfortran ./modulos/precisionMod.f90\
        ./modulos/mzranMod.f90\
        ./modulos/formatsMod.f90\
        ./modulos/constantsMod.f90\
        ./modulos/functionsMod.f90\
        ./modulos/parametersMod.f90\
        ./modulos/subroutinesMod.f90\
        ./modulos/potentialsMod.f90\
        ./modulos/observablesMod.f90\
        ./modulos/thermostatsMod.f90\
        ./modulos/updatePositionsMod.f90\
        ./modulos/parsingMod.f90\
        ./modulos/forcesMod.f90\
        ./modulos/writing2filesMod.f90\
        ./modulos/initializationsMod.f90\
        ./modulos/integratorsMod.f90\
        main.f90 -o run.exe -ffpe-trap=invalid,overflow,zero \
         -O3 -march=native -ftree-vectorize -ftree-vectorizer-verbose=2 \
         -O -Wall -fcheck=all -g -fbacktrace \
         -fopenmp -O2

echo
