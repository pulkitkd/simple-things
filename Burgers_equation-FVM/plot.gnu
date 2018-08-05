set term postscript
set out 'burgers.ps'
set key font ",16"
set xtics font ",20"
set ytics font ",20"
set xlabel font ",16"
set ylabel font ",16"

set xlabel 'x'
set ylabel 'u'
set xran[0:1]
set yran[-1.01:1.01]
plot 'burgersnum.dat', 'burgersexact.dat' with lines linecolor rgbcolor "red"
