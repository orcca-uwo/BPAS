#!/bin/bash

gnuplot << EOF
set term png
set output "TFT_experimental_data_4cores.png"

set size square 1,1
set title "TFT_experimental_data(4 cores)"


set xrange [13:24]
set yrange [-3:0.5]

set xtics 13,1,24
set ytics -3,0.2,0.5

set xlabel "size of input vector"
set ylabel "running time"
 
set grid

#set style line 

set key left top box

plot "tft_test_data.txt" using 1:2 title "tft python speedup" with linespoints pt 3 linecolor rgb "red","fft_test_data.txt" using 1:2 title "fft python speedup" with linespoints pt 3 linecolor rgb "blue","tft_parallel_test_data.txt" using 1:2 title "tft parallel python speedup" with linespoints pt 3 linecolor rgb "green"
EOF
