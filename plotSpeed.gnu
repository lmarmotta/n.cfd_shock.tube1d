#!/usr/bin/gnuplot

set term pdf monochrome size 15.0cm,9.0cm

set output "speed.pdf"

set grid

set xtics font "Times-Roman, 13"
set ytics font "Times-Roman, 13"

set xlabel "x position"
set ylabel "U Speed"
set title "U Speed plots - Analytical vs Numeric"

set offset graph 0.10,0.10,0.10,0.10

set offset graph 0.10,0.10,0.10,0.10

set border lw 2

set pointsize 0.3

set key font ",13"

plot "aspeed.out" u 1:2 with linespoints lt -1 lw 0.1 pt 4 title "Exact",\
     "velocityOutput.out" u 1:2 with lines lw 1.5 title "Numerical"
