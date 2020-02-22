#!/bin/bash

gnuplot << EOF
set term png
set output "TFT_experimental_data_maple.png"

set size square 1,1
set title "TFT_experimental_data"


set xrange [11:24]
set yrange [-3:3]

set xtics 11,1,24
set ytics -3,0.5,3

set xlabel "size of input vector"
set ylabel "running time"
 
#set grid

#set style line 

set key left top box

plot "maple_test_data.txt" using 1:2 title "tft maple speedup" with linespoints pt 3 linecolor rgb "green","tft_test_data.txt" using 1:2 title "tft python speedup" with linespoints pt 3 linecolor rgb "red"
EOF
