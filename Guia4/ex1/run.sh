./ex1.exe
(head -n 1 datos/temperature_functions.out && tail -n +2 datos/temperature_functions.out | sort -g) > datos/temperature_functions_sorted.out
gnuplot plot_observables.gp