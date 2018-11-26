set term postscript color enhanced 20
set output 'mtildemomentrates.ps'
#set colors classic

set xlabel 'Time (s)'
set ylabel 'Moment rate (*10^{18}Nm/s)'

#set yrange [0:1.6]

set title 'Mw6.2 24/08/2016'
plot 'mtildemomentrate0824.dat' u 1:($2/1.e18) t 'Kinematic inversion' w l lw 4 lc 2,\
'mtildemomentrate.dat' u 1:($2/1.e18) t 'Dynamic model' w l lw 4 lc 8
