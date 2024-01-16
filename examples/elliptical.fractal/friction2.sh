#!/bin/bash
#KORIGOVANA VERZE PRO NENULOVY DYNAMICKY KOEF. TRENI
rm slipcontour.dat nucleationzonecontour.dat
python >temp.out << END
from numpy import *
import matplotlib.pyplot as plt
f=open('inputfd3d.dat','r').readlines()
nxt,nzt=int(f[0].split()[0]),int(f[0].split()[2])
dh=float(f[1].split()[0])/1000.
f=open('input.dat','r').readlines()
L,W=float(f[17].split()[0])/1000.,float(f[17].split()[1])/1000.

dyn=0.4
coh=0.5

fricpar=loadtxt('result/friction.inp')
T0x=fricpar[:,0].reshape((nzt,nxt))/1.e6
T0z=fricpar[:,1].reshape((nzt,nxt))/1.e6
T0=T0x
Ts=fricpar[:,2].reshape((nzt,nxt))/1.e6
Dc=fricpar[:,3].reshape((nzt,nxt))
Mus=fricpar[:,4].reshape((nzt,nxt))
sd=loadtxt('result/stressdropX.res').reshape((nzt,nxt))/(-1.e6)
slip=loadtxt('result/slipX.res').reshape((nzt,nxt))
slipmax=slip.max()

x=linspace(dh/2.,L-dh/2.,nxt)
z=linspace(dh/2.,W-dh/2.,nzt)


f=open('slipcontour.dat','a')
for i in plt.contour(x,z,slip,[slipmax*.05]).collections[0].get_paths()[:]:
  savetxt(f,i.vertices)
  f.write('\n')
f.close()
f=open('nucleationzonecontour.dat','a')
for i in plt.contour(x,z,Ts+coh-T0,[0.]).collections[0].get_paths()[:]:
  savetxt(f,i.vertices)
  f.write('\n')
f.close()

#savetxt('slipcontour.dat',plt.contour(x,z,slip,[slipmax*.05]).collections[0].get_paths()[0].vertices)
#savetxt('nucleationzonecontour.dat',plt.contour(x,z,Ts+coh-T0,[0.]).collections[0].get_paths()[0].vertices)

f=open('friction.gnuplot.dat','w')
for j in range(nzt):
  f.write(' '.join([str(T0[j,i]-dyn*Ts[j,i]/Mus[j,i]) for i in range(nxt)])+'\n')
f.write('\n')
f.write('\n')
for j in range(nzt):
  f.write(' '.join([str(Mus[j,i]-dyn) for i in range(nxt)])+'\n')
f.write('\n')
f.write('\n')
for j in range(nzt):
  f.write(' '.join([str(Dc[j,i]) for i in range(nxt)])+'\n')
f.write('\n')
f.write('\n')
for j in range(nzt):
  f.write(' '.join([str(Ts[j,i]+coh-T0[j,i]) for i in range(nxt)])+'\n')
f.write('\n')
f.write('\n')
for j in range(nzt):
  f.write(' '.join([str(Dc[j,i]*(Ts[j,i]-dyn*Ts[j,i]/Mus[j,i])/2.) for i in range(nxt)])+'\n')
f.write('\n')
f.write('\n')
for j in range(nzt):
#  f.write(' '.join([str(min(Dc[j,i],slip[j,i])*Ts[j,i]/2.) for i in range(nxt)])+'\n')
  f.write(' '.join([str((Ts[j,i]-dyn*Ts[j,i]/Mus[j,i])/2.*(Dc[j,i]-max(0.,Dc[j,i]-slip[j,i])*(Dc[j,i]-slip[j,i])/Dc[j,i])) for i in range(nxt)])+'\n')
f.write('\n')
f.write('\n')
for j in range(nzt):
#  f.write(' '.join([str((slip[j,i]*Ts[j,i]-min(Dc[j,i],slip[j,i])*Ts[j,i])/2.) for i in range(nxt)])+'\n')
#  f.write(' '.join([str(Ts[j,i]/2.*(max(Dc[j,i],slip[j,i])-Dc[j,i])) for i in range(nxt)])+'\n')
  f.write(' '.join([str(((T0[j,i]-dyn*Ts[j,i]/Mus[j,i]+(Ts[j,i]-dyn*Ts[j,i]/Mus[j,i])*max(0.,1.-slip[j,i]/Dc[j,i]))*slip[j,i]-(Ts[j,i]-dyn*Ts[j,i]/Mus[j,i])*(Dc[j,i]-max(0.,Dc[j,i]-slip[j,i])**2/Dc[j,i]))/2.) for i in range(nxt)])+'\n')
f.write('\n')
f.write('\n')
for j in range(nzt):
  f.write(' '.join([str((Ts[j,i]-T0[j,i])/(T0[j,i]-dyn*Ts[j,i]/Mus[j,i])) for i in range(nxt)])+'\n')
f.write('\n')
f.write('\n')

rupt=loadtxt('result/ruptime.res').reshape((nzt,nxt))
f2=open('fractureengvsslip.dat','w')
for j in range(nzt):
  for i in range(nxt):
    f2.write(str((Ts[j,i]-dyn*Ts[j,i]/Mus[j,i])/2.*(Dc[j,i]-max(0.,Dc[j,i]-slip[j,i])*(Dc[j,i]-slip[j,i])/Dc[j,i]))+' '+str(rupt[j,i])+'\n')

print(' '.join([str(nxt),str(nzt),str(dh),str(L),str(W)]))
savetxt('inputmodel.meanvalues.dat',[mean(T0*slip)/mean(slip),mean(Ts*slip)/mean(slip),mean(Dc*slip)/mean(slip)])
END
read nxt nzt dh L W < temp.out

gnuplot << END
set term postscript color solid enhanced 12
set output 'frictionin2.ps'
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

set origin 0.,0.72
set title 'Initial stress (MPa)'
#set cbtics 10.
#set cbrange [0:7]
plot 'friction.gnuplot.dat' matrix u (\$1*$dh):(\$2*$dh):3 index 0 notitle w image,\
'slipcontour.dat' w l notitle lt -1 lw 1,\
'nucleationzonecontour.dat' w l notitle lt 3 lw 2,\
"aftershocks-onfault.dat" notitle w p ps .2 lc 9

set origin .5,.72
set title 'Static - Dynamic fric. coeff'
#set cbtics 0.2
#set cbrange [0:.25]
plot 'friction.gnuplot.dat' matrix u (\$1*$dh):(\$2*$dh):3 index 1 notitle w image,\
'slipcontour.dat' w l notitle lt -1 lw 1,\
'nucleationzonecontour.dat' w l notitle lt 3 lw 2,\
"aftershocks-onfault.dat" notitle w p ps .2 lc 9

set origin 0.,0.42
set title 'Dc (m)'
#set cbtics 0.05
#set cbrange [0:1]
plot 'friction.gnuplot.dat' matrix u (\$1*$dh):(\$2*$dh):3 index 2 notitle w image,\
'slipcontour.dat' w l notitle lt -1 lw 1,\
'nucleationzonecontour.dat' w l notitle lt 3 lw 2,\
"aftershocks-onfault.dat" notitle w p ps .2 lc 9
#set cbrange [:]

set origin 0.5,0.42
set title 'Strength excess (MPa)'
set cbtics 5.
set autoscale cb
set cbrange [0:20]
plot 'friction.gnuplot.dat' matrix u (\$1*$dh):(\$2*$dh):3 index 3 notitle w image,\
'slipcontour.dat' w l notitle lt -1 lw 1,\
'nucleationzonecontour.dat' w l notitle lt 3 lw 2,\
"aftershocks-onfault.dat" notitle w p ps .2 lc 9

set origin 0.,0.12
set title 'S-parameter'
set cbrange [0:5]
set cbtics 0.2
plot 'friction.gnuplot.dat' matrix u (\$1*$dh):(\$2*$dh):3 index 7 notitle w image,\
'slipcontour.dat' w l notitle lt -1 lw 1,\
'nucleationzonecontour.dat' w l notitle lt 3 lw 2,\
"aftershocks-onfault.dat" notitle w p ps .2 lc 9

set origin 0.5,0.12
#set title 'Breakdown work density (MJ/m^2)'
set title 'Fracture energy (MJ/m^2)'
set autoscale cb
set cbrange [0:2]
set cbtics 2.
plot 'friction.gnuplot.dat' matrix u (\$1*$dh):(\$2*$dh):3 index 4 notitle w image,\
'slipcontour.dat' w l notitle lt -1 lw 1,\
'nucleationzonecontour.dat' w l notitle lt 3 lw 2,\
"aftershocks-onfault.dat" notitle w p ps .2 lc 9

set origin 0.,-0.18
set title 'Dissipated breakdown work density (MJ/m^2)'
set autoscale cb
set cbrange [0:2]
set cbtics 2.
plot 'friction.gnuplot.dat' matrix u (\$1*$dh):(\$2*$dh):3 index 5 notitle w image,\
'slipcontour.dat' w l notitle lt -1 lw 1,\
'nucleationzonecontour.dat' w l notitle lt 3 lw 2,\
"aftershocks-onfault.dat" notitle w p ps .2 lc 9

set origin 0.5,-0.18
set title 'Radiated energy density (MJ/m^2)'
set cbrange [0:5]
#set autoscale cb
set cbtics 5.
plot 'friction.gnuplot.dat' matrix u (\$1*$dh):(\$2*$dh):3 index 6 notitle w image,\
'slipcontour.dat' w l notitle lt -1 lw 1,\
'nucleationzonecontour.dat' w l notitle lt 3 lw 2,\
"aftershocks-onfault.dat" notitle w p ps .2 lc 9

END
