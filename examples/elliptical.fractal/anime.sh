#!/bin/bash
FRAMES=560
REAL=1
mkdir tmp
rm -f tmp/*

NT=$FRAMES
DT=0.025
NL=200
NW=100
L=40.
W=20.
source load_gmt5
MAX=`minmax -C mtildeX.dat | awk '{print $2*0.8}'`
echo $MAX

python << END
from numpy import *
sliprate=loadtxt('mtildeX.dat')
print(sliprate.size,sliprate.max())
sliprate=sliprate.reshape(($NW,$NL,$NT))
print(sliprate.shape)
f=open('mtilde.gnuplot.dat','w')
for k in range($NT):
  for j in range($NW):
    f.write(' '.join([str(sliprate[j,i,k]) for i in range($NL)])+'\n')
  f.write('\n')
  f.write('\n')
END

count=1000
for ((i = 1; i <= $FRAMES; i = i + 1)); 
do
count=`expr $count + 1`
gnuplot << END
set term png size 1280,960
set output 'tmp/$count.png'
set pm3d map corners2color c3
set palette defined ( 0 "white", 2 "skyblue", 3 "light-green", 6 "yellow", 10 "light-red" )
set xtics out scale .2
set ytics out scale .2
set size ratio -1
DL=$L/($NL-1)
DW=$W/($NW-1)
T=0
set cbrange [0:$MAX]
set xrange [0:$L]
set yrange [0:$W]
set xlabel 'Along strike (km)'
set ylabel 'Up dip (km)'
set cblabel 'Slip rate (m/s)'
#set cbtics 0.1
splot 'mtilde.gnuplot.dat' matrix u (\$1*DL):(\$2*DW):3 index ($REAL-1)*$FRAMES+$i-1 notitle w pm3d
#"aftershocks-onfault.dat" notitle w p pt 1 ps .4 lc 0,\
#'epic.dat' u 1:2:(0) notitle w p pt 3 ps 2 lc 4
END
done

for i in `ls -1 tmp | grep png`;
do
convert tmp/$i tmp/$i.jpg;
done

mencoder mf://tmp/*.png.jpg -mf w=1280:h=960:fps=40. -o anime.avi -ovc lavc -lavcopts vcodec=mpeg4:vhq:mbd=1:vbitrate=600000
