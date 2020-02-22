#!/bin/bash

gnuplot << EOF
set term png
set output "TFT_relax_data_22to24_4cores.png"

set size square 1,1
set title "TFT_experimental_data_22to24(4 cores)"


set xrange [21.7:24.3]
set yrange [-0.15:0.8]

set xtics 21.7,0.3,24.3
set ytics -0.15,0.05,0.8

set xlabel "size of input vector"
set ylabel "running time"
 
#set grid

#set style line 

set key left top box

plot "tft_relax_22to24.txt" using 1:2 title "tft relax 22to24 python speedup" with linespoints pt 3 linecolor rgb "green"
EOF
