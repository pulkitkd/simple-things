# set terminal postscript eps enhanced color font 'Helvetica,16'
set terminal pngcairo size 1280,760 enhanced font 'Verdana,14'

set key right top Left title ''

set border ls 0

set xrange [0:800]

set yrange [-0.4:2.0]

# set grid ytics mytics  
# draw lines for each ytics and mytics

# set mytics 4 
set xlabel "density"

set ylabel "grid points"

set title "density evolution"




do for [t=0:4999] {
  outfile = sprintf('plots/density%03.0f.png',t)
  set output outfile 
  plot 'results/density_'.t.'.txt'  using 1:2 with lp title "density-".t lc 1 lw 2
}
