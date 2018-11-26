dh=0.1

set term postscript color landscape enhanced 11
set output 'processBayes.meansigma1.ps'
set multiplot
set colors classic
set size .45,.29
set palette defined ( 0 "white", 2 "skyblue", 3 "light-green", 6 "yellow", 10 "light-red" )
set xtics out scale .2
set ytics out scale .2
set size ratio -1
set xrange [0:30.]
set yrange [0:14.]
set xlabel 'Along strike (km)'
set ylabel 'Up dip (km)'

set origin 0.,0.82
set title 'S-parameter'
set cbrange [0:1]
plot 'processBayes.meansigma.dat' matrix u ($1*dh):($2*dh):($3<0?0:$3) index 8 notitle w image,\
'processBayes.slipcontour.dat' index 1 w l notitle lt -1 lw 1

set origin 0.5,0.82
set title '2xSigma S-parameter'
set cbrange [0:.5]
plot 'processBayes.meansigma.dat' matrix u ($1*dh):($2*dh):3 index 9 notitle w image,\
'processBayes.slipcontour.dat' index 1 w l notitle lt -1 lw 1

