#!/bin/bash
python >temp.out << END
from numpy import *
import matplotlib.pyplot as plt
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
Mus=fricpar[:,4].reshape((nzt,nxt))
sd=loadtxt('result/stressdrop.res').reshape((nzt,nxt))/(-1.e6)
slip=loadtxt('result/slip.res').reshape((nzt,nxt))
slipmax=slip.max()

x=linspace(dh/2.,L-dh/2.,nxt)
z=linspace(dh/2.,W-dh/2.,nzt)
savetxt('slipcontour.dat',plt.contour(x,z,slip,[slipmax*.05]).collections[0].get_paths()[0].vertices)
savetxt('nucleationzonecontour.dat',plt.contour(x,z,Ts-T0,[0.]).collections[0].get_paths()[0].vertices)

f=open('friction.gnuplot.dat','w')
for j in range(nzt):
  f.write(' '.join([str(T0[j,i]) for i in range(nxt)])+'\n')
f.write('\n')
f.write('\n')
for j in range(nzt):
  f.write(' '.join([str(Mus[j,i]) for i in range(nxt)])+'\n')
f.write('\n')
f.write('\n')
for j in range(nzt):
  f.write(' '.join([str(Dc[j,i]) for i in range(nxt)])+'\n')
f.write('\n')
f.write('\n')
for j in range(nzt):
  f.write(' '.join([str(Ts[j,i]-T0[j,i]) for i in range(nxt)])+'\n')
f.write('\n')
f.write('\n')
for j in range(nzt):
  f.write(' '.join([str(Dc[j,i]*Ts[j,i]/2.) for i in range(nxt)])+'\n')
f.write('\n')
f.write('\n')
for j in range(nzt):
#  f.write(' '.join([str(min(Dc[j,i],slip[j,i])*Ts[j,i]/2.) for i in range(nxt)])+'\n')
  f.write(' '.join([str(Ts[j,i]/2.*(Dc[j,i]-max(0.,Dc[j,i]-slip[j,i])*(Dc[j,i]-slip[j,i])/Dc[j,i])) for i in range(nxt)])+'\n')
f.write('\n')
f.write('\n')
for j in range(nzt):
#  f.write(' '.join([str((slip[j,i]*Ts[j,i]-min(Dc[j,i],slip[j,i])*Ts[j,i])/2.) for i in range(nxt)])+'\n')
#  f.write(' '.join([str(Ts[j,i]/2.*(max(Dc[j,i],slip[j,i])-Dc[j,i])) for i in range(nxt)])+'\n')
  f.write(' '.join([str(((T0[j,i]+Ts[j,i]*max(0.,1.-slip[j,i]/Dc[j,i]))*slip[j,i]-Ts[j,i]*(Dc[j,i]-max(0.,Dc[j,i]-slip[j,i])**2/Dc[j,i]))/2.) for i in range(nxt)])+'\n')
f.write('\n')
f.write('\n')
for j in range(nzt):
  f.write(' '.join([str((Ts[j,i]-T0[j,i])/T0[j,i]) for i in range(nxt)])+'\n')
f.write('\n')
f.write('\n')
print(' '.join([str(nxt),str(nzt),str(dh),str(L),str(W)]))
savetxt('inputmodel.meanvalues.dat',[mean(T0*slip)/mean(slip),mean(Ts*slip)/mean(slip),mean(Dc*slip)/mean(slip)])
END
read nxt nzt dh L W < temp.out

gnuplot << END
set term postscript color solid enhanced 12
set output 'frictionin2.ps'
set multiplot
set size 0.4,0.4
set palette defined ( 0 "white", 2 "skyblue", 3 "light-green", 6 "yellow", 10 "light-red" )
set xtics out scale .2
set ytics out scale .2
set size ratio -1
#set cbrange [0:$MAX]
set xrange [0:$nxt*$dh]
set yrange [0:$nzt*$dh]
set xlabel 'Along strike (km)'
set ylabel 'Up dip (km)'

set origin 0.,0.77
set title 'Prestress (MPa)'
#set cbtics 0.1
set cbrange [:]
plot 'friction.gnuplot.dat' matrix u (\$1*$dh):(\$2*$dh):3 index 0 notitle w image,\
'slipcontour.dat' w l notitle lt -1 lw 1,\
'nucleationzonecontour.dat' w l notitle lt 3 lw 2,\
"aftershocks-onfault.dat" notitle w p ps .2 lc 9

set origin .5,.77
set title 'Static - Dynamic fric. coeff'
#set cbtics 0.1
set cbrange [0:.6]
plot 'friction.gnuplot.dat' matrix u (\$1*$dh):(\$2*$dh):3 index 1 notitle w image,\
'slipcontour.dat' w l notitle lt -1 lw 1,\
'nucleationzonecontour.dat' w l notitle lt 3 lw 2,\
"aftershocks-onfault.dat" notitle w p ps .2 lc 9

set origin 0.,0.47
set title 'Dc (m)'
#set cbtics 0.1
set cbrange [.15:.8]
plot 'friction.gnuplot.dat' matrix u (\$1*$dh):(\$2*$dh):3 index 2 notitle w image,\
'slipcontour.dat' w l notitle lt -1 lw 1,\
'nucleationzonecontour.dat' w l notitle lt 3 lw 2,\
"aftershocks-onfault.dat" notitle w p ps .2 lc 9
set cbrange [:]

set origin 0.5,0.47
set title 'Strength excess (MPa)'
#set cbtics 0.1
set autoscale cb
set cbrange [-2:10]
plot 'friction.gnuplot.dat' matrix u (\$1*$dh):(\$2*$dh):(\$3) index 3 notitle w image,\
'slipcontour.dat' w l notitle lt -1 lw 1,\
'nucleationzonecontour.dat' w l notitle lt 3 lw 2,\
"aftershocks-onfault.dat" notitle w p ps .2 lc 9

set origin 0.,0.17
set title 'S-parameter'
set cbrange [0:1]
#set cbtics 0.1
plot 'friction.gnuplot.dat' matrix u (\$1*$dh):(\$2*$dh):3 index 7 notitle w image,\
'slipcontour.dat' w l notitle lt -1 lw 1,\
'nucleationzonecontour.dat' w l notitle lt 3 lw 2,\
"aftershocks-onfault.dat" notitle w p ps .2 lc 9

set origin 0.5,0.17
set title 'Breakdown work density (MJ/m^2)'
set autoscale cb
set cbrange [0:3]
#set cbtics 0.1
plot 'friction.gnuplot.dat' matrix u (\$1*$dh):(\$2*$dh):3 index 4 notitle w image,\
'slipcontour.dat' w l notitle lt -1 lw 1,\
'nucleationzonecontour.dat' w l notitle lt 3 lw 2,\
"aftershocks-onfault.dat" notitle w p ps .2 lc 9

set origin 0.,-0.13
set title 'Dissipated breakdown work density (MJ/m^2)'
set autoscale cb
set cbrange [0:3]
#set cbtics 0.1
plot 'friction.gnuplot.dat' matrix u (\$1*$dh):(\$2*$dh):3 index 5 notitle w image,\
'slipcontour.dat' w l notitle lt -1 lw 1,\
'nucleationzonecontour.dat' w l notitle lt 3 lw 2,\
"aftershocks-onfault.dat" notitle w p ps .2 lc 9

set origin 0.5,-0.13
set title 'Radiated energy density (MJ/m^2)'
set cbrange [0:6]
set autoscale cb
#set cbtics 0.1
plot 'friction.gnuplot.dat' matrix u (\$1*$dh):(\$2*$dh):3 index 6 notitle w image,\
'slipcontour.dat' w l notitle lt -1 lw 1,\
'nucleationzonecontour.dat' w l notitle lt 3 lw 2,\
"aftershocks-onfault.dat" notitle w p ps .2 lc 9

END
