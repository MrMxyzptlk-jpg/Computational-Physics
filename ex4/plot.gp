set terminal x11 

set title "Error relativo vs Numero de puntos"
set ylabel "Error Relativo"
set xlabel "N"
set grid
set logscale
set key box opaque 
set key top right

plot "./datos/datos_ex4.out" using 2:5 with linespoints title "Trapecio",\
    "./datos/datos_ex4.out" using 2:7 with linespoints title "Simpson",\
    "./datos/datos_ex4.out" using 2:9 with linespoints title "Gauss",\

set terminal pngcairo
set output "grafico.png"
replot

set terminal x11