set terminal png font "sans-serif,10"
set terminal qt
set key top center
plot "solAnalitica.dat" using 1:2 with lines, \
"solNumerica.dat" using 1:2 with points ps 1 pt 7
replot
