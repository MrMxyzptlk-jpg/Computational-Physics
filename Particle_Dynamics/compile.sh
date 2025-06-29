gfortran src/modulos/precisionMod.f90\
        src/modulos/mzranMod.f90\
        src/modulos/mzran_threadsafeMod.f90\
        src/modulos/randomMod.f90\
        src/modulos/formatsMod.f90\
        src/modulos/constantsMod.f90\
        src/modulos/parametersMod.f90\
        src/modulos/subroutinesMod.f90\
        src/modulos/potentialsMod.f90\
        src/modulos/observablesMod.f90\
        src/modulos/thermostatsMod.f90\
        src/modulos/updatePositionsMod.f90\
        src/modulos/parsingMod.f90\
        src/modulos/forcesMod.f90\
        src/modulos/writing2filesMod.f90\
        src/modulos/initializationsMod.f90\
        src/modulos/integratorsMod.f90\
        src/main.f90 -o run.exe -ffpe-trap=invalid,overflow,zero \
         -O3 -march=native -ftree-vectorize -ftree-vectorizer-verbose=2 \
         -O -Wall -fcheck=all -g -fbacktrace \
         -fopenmp -O2 \
         -Iexternal/FoX/objs/finclude -Lexternal/FoX/objs/lib -lFoX_wxml   # Include FoX library for parsing XML

echo
