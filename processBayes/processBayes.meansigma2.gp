dh=0.1

set term postscript color landscape enhanced 11
set output 'processBayes.meansigma2.ps'
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
set title 'Mean Slip (m)'
plot 'processBayes.meansigma2.dat' matrix u ($1*dh):($2*dh):3 index 0 notitle w image,\
'processBayes.slipcontour.dat' index 1 w l notitle lt -1 lw 1

set origin 0.5,0.82
set title '2xSigma Slip (m)'
#set cbrange [0:.1]
plot 'processBayes.meansigma2.dat' matrix u ($1*dh):($2*dh):3 index 1 notitle w image,\
'processBayes.slipcontour.dat' index 1 w l notitle lt -1 lw 1

set origin 0.,0.52
set title 'Mean Stress drop (MPa)'
plot 'processBayes.meansigma2.dat' matrix u ($1*dh):($2*dh):(-$3/1.e6) index 2 notitle w image,\
'processBayes.slipcontour.dat' index 1 w l notitle lt -1 lw 1

set origin 0.5,0.52
set title '2xSigma Stress drop (MPa)'
#set cbrange [0:.1]
plot 'processBayes.meansigma2.dat' matrix u ($1*dh):($2*dh):($3/1.e6) index 3 notitle w image,\
'processBayes.slipcontour.dat' index 1 w l notitle lt -1 lw 1

set origin 0.,0.22
set title 'Mean Rise time (s)'
plot 'processBayes.meansigma2.dat' matrix u ($1*dh):($2*dh):3 index 4 notitle w image,\
'processBayes.slipcontour.dat' index 1 w l notitle lt -1 lw 1

set origin 0.5,0.22
set title '2xSigma Rise time (s)'
#set cbrange [0:.1]
plot 'processBayes.meansigma2.dat' matrix u ($1*dh):($2*dh):3 index 5 notitle w image,\
'processBayes.slipcontour.dat' index 1 w l notitle lt -1 lw 1

set origin 0.,-0.08
set title 'Mean Rupture time (s)'
plot 'processBayes.meansigma2.dat' matrix u ($1*dh):($2*dh):3 index 6 notitle w image,\
'processBayes.slipcontour.dat' index 1 w l notitle lt -1 lw 1

set origin 0.5,-0.08
set title '2xSigma Rupture time (s)'
#set cbrange [0:.1]
plot 'processBayes.meansigma2.dat' matrix u ($1*dh):($2*dh):3 index 7 notitle w image,\
'processBayes.slipcontour.dat' index 1 w l notitle lt -1 lw 1
