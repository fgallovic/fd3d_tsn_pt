#!/bin/bash
python >temp.out << END
from numpy import *
f=open('inputfd3d.dat','r').readlines()
nxt,nzt=int(f[0].split()[0]),int(f[0].split()[2])
dh=float(f[1].split()[0])/1000.
f=open('input.dat','r').readlines()
L,W=float(f[17].split()[0])/1000.,float(f[17].split()[1])/1000.

fricpar=loadtxt('result/friction.inp')
T0x=fricpar[:,0].reshape((nzt,nxt))/1.e6
T0z=fricpar[:,1].reshape((nzt,nxt))/1.e6
T0=T0x
Ts=fricpar[:,2].reshape((nzt,nxt))/1.e6
Dc=fricpar[:,3].reshape((nzt,nxt))
f=open('friction.gnuplot.dat','w')
for j in range(nzt):
  f.write(' '.join([str(T0[j,i]) for i in range(nxt)])+'\n')
f.write('\n')
f.write('\n')
for j in range(nzt):
  f.write(' '.join([str(Ts[j,i]) for i in range(nxt)])+'\n')
f.write('\n')
f.write('\n')
for j in range(nzt):
  f.write(' '.join([str(Dc[j,i]) for i in range(nxt)])+'\n')
f.write('\n')
f.write('\n')
mu1=loadtxt('result/vmodel.inp')
mu=mu1.reshape((nzt,nxt))/1.e9
for j in range(nzt):
  f.write(' '.join([str(mu[j,i]) for i in range(nxt)])+'\n')
print(' '.join([str(nxt),str(nzt),str(dh),str(L),str(W)]))
END
read nxt nzt dh L W < temp.out

gnuplot << END
set term postscript color 12
set output 'friction1.ps'
set multiplot
set size 0.3,0.3
set palette defined ( 0 "white", 2 "skyblue", 3 "light-green", 6 "yellow", 10 "light-red" )
set xtics out scale .2
set ytics out scale .2
set size ratio -1
#set cbrange [0:$MAX]
set xrange [0:$nxt*$dh]
set yrange [0:$nzt*$dh]
set xlabel 'Along strike (km)'
set ylabel 'Up dip (km)'

set origin 0.,0.5
set title 'Prestress (MPa)'
#set cbtics 0.1
#set cbrange [-2:9]
plot 'friction.gnuplot.dat' matrix u (\$1*$dh):(\$2*$dh):3 index 0 notitle w image,\
"aftershocks-onfault.dat" notitle w p ps .5 lc 10

set origin .5,.5
set title 'Strength (MPa)'
#set cbtics 0.1
#set cbrange [0:2.5]
plot 'friction.gnuplot.dat' matrix u (\$1*$dh):(\$2*$dh):3 index 1 notitle w image,\
"aftershocks-onfault.dat" notitle w p ps .5 lc 10

set origin 0.,0.2
set title 'Dc (m)'
#set cbtics 0.1
plot 'friction.gnuplot.dat' matrix u (\$1*$dh):(\$2*$dh):3 index 2 notitle w image,\
"aftershocks-onfault.dat" notitle w p ps .5 lc 10

set origin 0.5,0.2
set title 'Mu (GPa)'
#set cbtics 0.1
plot 'friction.gnuplot.dat' matrix u (\$1*$dh):(\$2*$dh):3 index 3 notitle w image,\
"aftershocks-onfault.dat" notitle w p ps .5 lc 10

END
