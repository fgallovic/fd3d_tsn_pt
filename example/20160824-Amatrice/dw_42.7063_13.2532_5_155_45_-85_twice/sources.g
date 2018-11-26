set term postscript color
set output 'sources.ps'
set size ratio -1
#set xrange [-50:50]
#set yrange [-50:50]
plot 'sources.dat' u 3:2 w p,'fault.dat' u 2:1 w l,\
'stations.dat' u 2:1  w p,'stations.dat' u 2:1:4 w labels
