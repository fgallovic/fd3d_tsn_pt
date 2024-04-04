set term postscript color solid enhanced 16
set output 'processBayes.MomentRates.ps'
#set size .5,.5

set palette grey negative
set cbrange [0:1]
set xlabel 'Time (s)'
set ylabel 'Moment rate (Nm)'
set cblabel 'Normalized PDF'
plot 'processBayes.MomentRates.dat' w l notitle lc palette lw 5
#plot 'processBayes.MomentRates.dat' u 1:2:($3>.99?$3:1/0) w l notitle lc palette lt 9 lw 5
