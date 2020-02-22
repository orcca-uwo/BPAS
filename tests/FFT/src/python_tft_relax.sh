#!/bin/bash

gnuplot << EOF
set term png
set output "TFT_relax_data_4cores.png"

set size square 1,1
set title "TFT_experimental_data(4 cores)"


set xrange [15:17]
set yrange [-3:-1.0]

set xtics 15,0.3,17
set ytics -3,0.1,-1.0

set xlabel "size of input vector"
set ylabel "running time"
 
#set grid

#set style line 

set key left top box

plot "python_tft_relax.txt" using 1:2 title "tft relax python speedup" with linespoints pt 3 linecolor rgb "red"
EOF
