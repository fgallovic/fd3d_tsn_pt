#!/bin/bash
#nxt=180
#nzt=60
#dh=0.2
#L=36.
#W=12.
#ntfd=1600
#dt=0.005
#MAX=`minmax -C mtilde.dat | awk '{print $2*0.8}'`
echo $MAX

python >temp.out << END
from numpy import *
import matplotlib.pyplot as plt
f=open('inputfd3d.dat','r').readlines()
nxt,nzt=int(f[0].split()[0]),int(f[0].split()[2])
dh=float(f[1].split()[0])/1000.
f=open('input.dat','r').readlines()
L,W=float(f[17].split()[0])/1000.,float(f[17].split()[1])/1000.

f=open('outputmodel.gnuplot.dat','w')
slip=loadtxt('result/slip.res').reshape((nzt,nxt))
sd=loadtxt('result/stressdrop.res').reshape((nzt,nxt))/(-1.e6)
rupt=loadtxt('result/ruptime.res').reshape((nzt,nxt))
rist=loadtxt('result/risetime.res').reshape((nzt,nxt))
Dc=loadtxt('result/friction.inp')[:,2].reshape((nzt,nxt))
slipmax=slip.max()

x=linspace(dh/2.,L-dh/2.,nxt)
z=linspace(dh/2.,W-dh/2.,nzt)
savetxt('slipcontour.dat',plt.contour(x,z,slip,[slipmax*.05]).collections[0].get_paths()[0].vertices)

for j in range(nzt):
  f.write(' '.join([str(slip[j,i]) for i in range(nxt)])+'\n')
f.write('\n')
f.write('\n')
for j in range(nzt):
  f.write(' '.join([str(sd[j,i]) for i in range(nxt)])+'\n')
f.write('\n')
f.write('\n')
for j in range(nzt):
  f.write(' '.join([str(rupt[j,i]) for i in range(nxt)])+'\n')
f.write('\n')
f.write('\n')
for j in range(nzt):
  f.write(' '.join([str(rist[j,i]) for i in range(nxt)])+'\n')
print(' '.join([str(nxt),str(nzt),str(dh),str(L),str(W)]))
END
read nxt nzt dh L W < temp.out

gnuplot << END
set term postscript color 12
set output 'outputmodel.ps'
set multiplot
set size 0.5,0.5
set palette defined ( 0 "white", 2 "skyblue", 3 "light-green", 6 "yellow", 10 "light-red" )
set xtics out scale .2
set ytics out scale .2
set size ratio -1
set cbrange [0:]
set xrange [0:$nxt*$dh]
set yrange [0:$nzt*$dh]
set xlabel 'Along strike (km)'
set ylabel 'Up dip (km)'

set origin 0.,0.5
set title 'Slip (m)'
#set cbtics 0.1
plot 'outputmodel.gnuplot.dat' matrix u (\$1*$dh):(\$2*$dh):3 index 0 notitle w image,\
'slipcontour.dat' w l notitle lt -1 lw 1,\
"aftershocks-onfault.dat" notitle w p ps .1 lc 10

set origin .5,.5
set title 'Stress drop (MPa)'
#set cbtics 0.1
plot 'outputmodel.gnuplot.dat' matrix u (\$1*$dh):(\$2*$dh):3 index 1 notitle w image,\
'slipcontour.dat' w l notitle lt -1 lw 1,\
"aftershocks-onfault.dat" notitle w p ps .1 lc 10

set origin 0.,0.2
set title 'Rupture time (s)'
#set cbtics 0.1
plot 'outputmodel.gnuplot.dat' matrix u (\$1*$dh):(\$2*$dh):3 index 2 notitle w image,\
'slipcontour.dat' w l notitle lt -1 lw 1,\
"aftershocks-onfault.dat" notitle w p ps .1 lc 10

set origin 0.5,0.2
set title 'Risetime (s)'
#set cbtics 0.1
plot 'outputmodel.gnuplot.dat' matrix u (\$1*$dh):(\$2*$dh):3 index 3 notitle w image,\
'slipcontour.dat' w l notitle lt -1 lw 1,\
"aftershocks-onfault.dat" notitle w p ps .1 lc 10

END
