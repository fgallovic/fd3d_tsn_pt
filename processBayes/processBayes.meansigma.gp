dh=0.1

set term postscript color landscape enhanced 11
set output 'processBayes.meansigma.ps'
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
set title 'Mean Prestress (MPa)'
set cbrange [0:12]
plot 'processBayes.meansigma.dat' matrix u ($1*dh):($2*dh):3 index 0 notitle w image,\
'processBayes.slipcontour.dat' index 1 w l notitle lt -1 lw 1

set origin 0.5,0.82
set title '2xSigma Prestress (MPa)'
set cbrange [0:3]
plot 'processBayes.meansigma.dat' matrix u ($1*dh):($2*dh):3 index 1 notitle w image,\
'processBayes.slipcontour.dat' index 1 w l notitle lt -1 lw 1

set origin 0.,0.52
set title 'Mean Static - Dynamic fric. coeff.'
set cbrange [0:.6]
plot 'processBayes.meansigma.dat' matrix u ($1*dh):($2*dh):3 index 2 notitle w image,\
'processBayes.slipcontour.dat' index 1 w l notitle lt -1 lw 1

set origin 0.5,0.52
set title '2xSigma Static - Dynamic fric. coeff.'
set cbrange [0:.1]
plot 'processBayes.meansigma.dat' matrix u ($1*dh):($2*dh):3 index 3 notitle w image,\
'processBayes.slipcontour.dat' index 1 w l notitle lt -1 lw 1

set origin 0.,0.22
set title 'Mean Strength excess (MPa)'
set cbrange [-2:10]
plot 'processBayes.meansigma.dat' matrix u ($1*dh):($2*dh):3 index 4 notitle w image,\
'processBayes.slipcontour.dat' index 1 w l notitle lt -1 lw 1

set origin 0.5,0.22
set title '2xSigma Strength excess (MPa)'
set cbrange [0:4.]
plot 'processBayes.meansigma.dat' matrix u ($1*dh):($2*dh):3 index 5 notitle w image,\
'processBayes.slipcontour.dat' index 1 w l notitle lt -1 lw 1

set origin 0.,-0.08
set title 'Mean Dc (m)'
set autoscale cb
set cbrange [0.15:0.8]
plot 'processBayes.meansigma.dat' matrix u ($1*dh):($2*dh):3 index 6 notitle w image,\
'processBayes.slipcontour.dat' index 1 w l notitle lt -1 lw 1

set origin 0.5,-0.08
set title '2xSigma Dc (m)'
set cbrange [0:.15]
plot 'processBayes.meansigma.dat' matrix u ($1*dh):($2*dh):3 index 7 notitle w image,\
'processBayes.slipcontour.dat' index 1 w l notitle lt -1 lw 1
