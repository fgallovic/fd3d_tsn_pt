#!/bin/bash
rm slipcontour.dat nucleationzonecontour.dat
python << END >temp.out 
from numpy import *
import matplotlib.pyplot as plt
f=open('inputfd3d.dat','r').readlines()
nxt,nzt=int(f[0].split()[0]),int(f[0].split()[2])
dh=float(f[1].split()[0])/1000.
f=open('input.dat','r').readlines()
L,W=float(f[17].split()[0])/1000.,float(f[17].split()[1])/1000.

f=open('outputmodel.gnuplot.dat','w')
slip=loadtxt('result/slipX.res').reshape((nzt,nxt))
sd=loadtxt('result/stressdropX.res').reshape((nzt,nxt))/(-1.e6)
rupt=loadtxt('result/ruptime.res').reshape((nzt,nxt))
ruptvel=loadtxt('result/rvel.res').reshape((nzt,nxt))
rist=loadtxt('result/risetime.res').reshape((nzt,nxt))
fricpar=loadtxt('result/friction.inp')
T0x=fricpar[:,0].reshape((nzt,nxt))/1.e6
T0z=fricpar[:,1].reshape((nzt,nxt))/1.e6
T0=T0x
Ts=fricpar[:,2].reshape((nzt,nxt))/1.e6
slipmax=slip.max()
slipmean=slip.sum()/(slip>0.001).sum()
coh=0.5

x=linspace(dh/2.,L-dh/2.,nxt)
z=linspace(dh/2.,W-dh/2.,nzt)

f1=open('slipcontour.dat','a')
for i in plt.contour(x,z,slip,[slipmax*.05]).collections[0].get_paths()[:]:
  savetxt(f1,i.vertices)
  f1.write('\n')
f1.close()
f1=open('nucleationzonecontour.dat','a')
for i in plt.contour(x,z,Ts+coh-T0,[0.]).collections[0].get_paths()[:]:
  savetxt(f1,i.vertices)
  f1.write('\n')
f1.close()

#savetxt('slipcontour.dat',plt.contour(x,z,slip,[slipmax*.05]).collections[0].get_paths()[0].vertices)
#savetxt('nucleationzonecontour.dat',plt.contour(x,z,Ts+coh-T0,[0.]).collections[0].get_paths()[0].vertices)

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
f.write('\n')
f.write('\n')
for j in range(nzt):
  f.write(' '.join([str(ruptvel[j,i]/1000.) for i in range(nxt)])+'\n')

print(' '.join([str(nxt),str(nzt),str(dh),str(L),str(W),str(slipmean)]))
END
read nxt nzt dh L W < temp.out

gnuplot << END
set term postscript color solid enhanced 12
set output 'outputmodel.ps'
set multiplot
set size 0.35,0.35
set palette defined ( 0 "white", 2 "skyblue", 3 "light-green", 6 "yellow", 10 "light-red" )
set xtics out scale .2
set ytics out scale .2
set size ratio -1
set xrange [0:$nxt*$dh]
set yrange [0:$nzt*$dh]
set xlabel 'Along strike (km)'
set ylabel 'Up dip (km)'

set origin 0.,0.55
set title 'Stress drop (MPa)'
set cbtics 10.
set cbrange [0:10]
plot 'outputmodel.gnuplot.dat' matrix u (\$1*$dh):(\$2*$dh):(\$3) index 1 notitle w image,\
'slipcontour.dat' w l notitle lt -1 lw 1,\
'nucleationzonecontour.dat' w l notitle lt 3 lw 2,\
"aftershocks-onfault.dat" notitle w p ps .1 lc 9


set origin .5,.55
set title 'Risetime (s)'
set cbrange [0:15]
#set cbtics 0.1
plot 'outputmodel.gnuplot.dat' matrix u (\$1*$dh):(\$2*$dh):3 index 3 notitle w image,\
'slipcontour.dat' w l notitle lt -1 lw 1,\
'nucleationzonecontour.dat' w l notitle lt 3 lw 2,\
"aftershocks-onfault.dat" notitle w p ps .1 lc 9

set origin 0.,0.2
set title 'Rupture time (s)'
set cbtics 5.
set cbrange [0:15]
plot 'outputmodel.gnuplot.dat' matrix u (\$1*$dh):(\$2*$dh):3 index 2 notitle w image,\
'slipcontour.dat' w l notitle lt -1 lw 1,\
'nucleationzonecontour.dat' w l notitle lt 3 lw 2,\
"aftershocks-onfault.dat" notitle w p ps .2 lc 9

set origin 0.5,0.2
set title 'Rupture speed (km/s)'
set cbtics 1.
set cbrange [0:3.5]
#plot 'outputmodel.gnuplot.dat' matrix u (\$1*$dh):(\$2*$dh):(\$3<5?\$3:(5)) index 4 notitle w image
plot 'outputmodel.gnuplot.dat' matrix u (\$1*$dh):(\$2*$dh):3 index 4 notitle w image,\
'slipcontour.dat' w l notitle lt -1 lw 1,\
'nucleationzonecontour.dat' w l notitle lt 3 lw 2,\
"aftershocks-onfault.dat" notitle w p ps .1 lc 9
set autoscale cb

set origin 0.,-0.15
set title 'Slip (m)'
set cbtics 1
set cbrange [0:5.]
plot 'outputmodel.gnuplot.dat' matrix u (\$1*$dh):(\$2*$dh):3 index 0 notitle w image,\
'slipcontour.dat' w l notitle lt -1 lw 1,\
'nucleationzonecontour.dat' w l notitle lt 3 lw 2,\
"aftershocks-onfault.dat" notitle w p ps .1 lc 9

END
