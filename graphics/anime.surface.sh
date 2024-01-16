#!/bin/bash
NT=800
DT=0.025
NL=200
NW=400
L=20.
W=40.
FRAMES=$NT

mkdir tmp
rm -f tmp/*
python >temp.out << END
import numpy as np
S=np.loadtxt('result/seisoutV.surface.gnuplot.dat')
for k in range($NT):
  np.savetxt('tmp/'+str(1000+k+1)+'.dat',S[k*$NW:(k+1)*$NW])
print(abs(S).max()*0.8)
END
read MAX < temp.out
echo $MAX

count=1000
for ((i = 1; i <= $FRAMES; i = i + 1)); 
do
count=`expr $count + 1`
gnuplot << END
set term png size 1280,960
set output 'tmp/$count.png'
set palette defined ( -1 "blue", 0 "white", 1 "red" )
set xtics out scale .2
set ytics out scale .2
set size ratio -1
DL=$L/($NL)
DW=$W/($NW)
T=0
set cbrange [-$MAX:$MAX]
set xrange [0:$L]
set yrange [0:$W]
set xlabel 'North (km)'
set ylabel 'East (km)'
set cblabel 'Velocity (m/s)'
#set cbtics 0.1
plot 'tmp/$count.dat' matrix u (\$1*DL):(\$2*DW):3 notitle w image
#'epic.dat' u 1:2:(0) notitle w p pt 3 ps 2 lc 4
END
done

for i in `ls -1 tmp | grep png`;
do
convert tmp/$i tmp/$i.jpg;
done

mencoder mf://tmp/*.png.jpg -mf w=1280:h=960:fps=40. -o anime.surface.avi -ovc lavc -lavcopts vcodec=mpeg4:vhq:mbd=1:vbitrate=600000
