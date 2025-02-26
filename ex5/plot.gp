set terminal x11 
set xlabel "t"
set ylabel "y(t)"
#set xrange [-0.01:1.01]
set grid
set key box opaque 
set key top right

set terminal pngcairo
set output "grafico1.png"
set title "funcion real y aproximaciones con h=0.1"
file_in = "./datos/datos_ex5_1.out"
           
plot file_in using 1:3 with linespoints title "Euler", file_in using 1:5 with linespoints title "RK2",  \
    file_in using 1:7 with linespoints title "RK4", \
    exp(-x**2/2) title "primitiva"


set terminal pngcairo
set output "grafico2.png"

set title "funcion real y aproximaciones con h=0.01"
file_in = "./datos/datos_ex5_2.out"
           
plot file_in using 1:3 with linespoints title "Euler", file_in using 1:5 with linespoints title "RK2",  \
    file_in using 1:7 with linespoints title "RK4", \
    exp(-x**2/2) title "primitiva"


set terminal pngcairo
set output "grafico3.png"

set title "funcion real y aproximaciones con h=0.001"
file_in = "./datos/datos_ex5_3.out"
           
plot file_in using 1:3 with lines title "Euler", file_in using 1:5 with lines title "RK2",  \
    file_in using 1:7 with lines title "RK4", \
    exp(-x**2/2) title "primitiva"


set terminal pngcairo
set output "grafico4.png"

set title "funcion real y aproximaciones con h=0.0001"
file_in = "./datos/datos_ex5_4.out"
           
plot file_in using 1:3 with lines title "Euler", file_in using 1:5 with lines title "RK2",  \
    file_in using 1:7 with lines title "RK4", \
    exp(-x**2/2) title "primitiva"


set terminal pngcairo
set output "grafico5.png"

set title "funcion real y aproximaciones con h=0.00001"
file_in = "./datos/datos_ex5_5.out"
           
plot file_in using 1:3 with lines title "Euler", file_in using 1:5 with lines title "RK2",  \
    file_in using 1:7 with lines title "RK4", \
    exp(-x**2/2) title "primitiva"


set terminal x11 

set terminal pngcairo
set output "times.png"
set xlabel "tiempo [s]"
set ylabel "Precisión"
set logscale x
set title "Tiempos de cálculo en función de la precisión alcanzada"
file_in = "./datos/times.out"

plot file_in using 2:1 with lines title "Euler", file_in using 4:3 with lines title "RK2",  \
    file_in using 6:5 with lines title "RK4"
set terminal x11 

