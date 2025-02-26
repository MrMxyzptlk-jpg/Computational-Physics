set terminal x11

set title "Error versus n"
set xlabel "n"
set ylabel "Error absoluto "
set grid
set logscale
set key box opaque 
#set key top center
           

plot "./datos/datos_ex3.out" using 2:5 with linespoints title "Trapecio","./datos/datos_ex3.out" using 2:7 with \
    linespoints title "simpson"
    
set terminal pngcairo
set output "./grafico.png"
replot

set terminal x11