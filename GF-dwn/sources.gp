set term postscript 6
set output 'sources.ps'
set size ratio -1
#set xrange [-30:30]
#set yrange [-30:30]
plot 'sources.dat' u 3:2 w p,\
'fault.dat' u ($2/1000.):($1/1000.) notitle w l,\
'stations.dat' u 2:1 notitle,\
'stations.dat' u 2:1:4 notitle w labels
