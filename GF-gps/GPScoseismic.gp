 set term postscript portrait color solid enh
 set output "gpscoseismic.ps"
set multiplot
set size .5,1
set size ratio -1
norm=1000.
set xrange [-50:50]
set yrange [-50:50]
set palette defined ( 0 "white", 2 "skyblue", 3 "light-green", 6 "yellow", 10 "light-red" )
#set cbrange [0:1]

set origin 0,0
plot 'sources.dat' u 3:2 w p ps .05,\
'stations-GPS.dat' u 2:1:($4*norm):($3*norm) t 'Data GPS' w vectors lt -1 lw 2,\
'mtildeslip2D-sGPS.out' u 2:1:($4*norm):($3*norm) notitle w vectors lw 1,\
'stations-GPS.dat' u 2:1:6 w labels

set origin .5,0
plot 'sources.dat' u 3:2 w p ps .05,\
'stations-GPS.dat' u 2:1:(0):($5*norm) t 'Data GPS' w vectors lt -1 lw 2,\
'mtildeslip2D-sGPS.out' u 2:1:(0):($5*norm) notitle w vectors lw 1,\
'stations-GPS.dat' u 2:1:6 w labels

unset multiplot
