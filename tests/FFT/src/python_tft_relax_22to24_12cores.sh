#!/bin/bash

gnuplot << EOF
set term png
set output "TFT_relax_data_22to24_12cores.png"

set size square 1,1
set title "TFT_experimental_data_22to24(12 cores)"


set xrange [21.7:24.3]
set yrange [-1.15:-0.3]

set xtics 21.7,0.3,24.3
set ytics -1.15,0.1,-0.3

set xlabel "size of input vector"
set ylabel "running time"
 
#set grid

#set style line 

set key left top box

plot "tft_relax_22to24_12cores.txt" using 1:2 title "tft relax 22to24 python speedup" with linespoints pt 3 linecolor rgb "green"
EOF
