set terminal x11

set title "Error versus h"
set xlabel "incremento h"
set ylabel "Error de df(x)/dx"
set grid
set logscale 

plot "./datos/datos_ex1.out" using 1:5 with linespoints title ""
set terminal pngcairo
set output "grafico.png"
replot

set terminal x11