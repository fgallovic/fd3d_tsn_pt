#!/bin/bash
nxt=200
nzt=200
dh=0.01
L=2.
W=2.

python << END
from numpy import *
import matplotlib.pyplot as plt
slip=loadtxt('processBayes.slipmodels.dat')
misfits=loadtxt('processBayes.dat')[:,0]
maxmisfit=min(misfits)
NM=int(slip.size/$nxt/$nzt)
sortmisfit=argsort(-misfits)
slip=slip.reshape(($nzt,$nxt,NM))
meanslip=slip.mean(axis=2)
stdslip=2*slip.std(axis=2)
stdslipp=meanslip+stdslip
stdslipm=meanslip-stdslip

x=linspace($dh/2.,$L-$dh/2.,$nxt)
z=linspace($dh/2.,$W-$dh/2.,$nzt)
f= open('processBayes.slipcontour.dat','w')
for i in sortmisfit:
  slip1=slip[:,:,i]
  slipmax=slip1.max()
  cnt=plt.contour(x,z,slip1,[slipmax*.05]).collections[0].get_paths()[0].vertices
  f.write('\n'.join([str(s[0])+' '+str(s[1])+' '+str(exp(-(misfits[i]-maxmisfit))) for s in cnt]))
  f.write('\n\n')
f.write('\n\n')

for i in plt.contour(x,z,meanslip,[meanslip.max()*.05]).collections[0].get_paths()[:]:
  savetxt(f,i.vertices)
  f.write('\n')
f.write('\n\n')
for i in plt.contour(x,z,stdslipp,[stdslipp.max()*.05]).collections[0].get_paths()[:]:
  savetxt(f,i.vertices)
  f.write('\n')
f.write('\n\n')
for i in plt.contour(x,z,stdslipm,[stdslipm.max()*.05]).collections[0].get_paths()[:]:
  savetxt(f,i.vertices)
  f.write('\n')
f.write('\n\n')

END

gnuplot << END
set term postscript color solid enhanced 8
set output 'processBayes.slipmodels.ps'
set multiplot
set size 0.5,0.5
set xtics out scale .2
set ytics out scale .2
set size ratio -1
set xrange [0:$nxt*$dh]
set yrange [0:$nzt*$dh]
set xlabel 'Along strike (km)'
set ylabel 'Up dip (km)'
set cblabel 'Normalized PDF'
set cbrange [0:]

set origin 0.,0.5
set title 'Slip contours'

set palette grey negative
plot 'processBayes.slipcontour.dat' index 0 w l notitle lt 9 lw .1 lc palette,\
'processBayes.slipcontour.dat' index 1 w l notitle lt 9 lw 1,\
'processBayes.slipcontour.dat' index 2 w l notitle lt 9 lw .5,\
'processBayes.slipcontour.dat' index 3 w l notitle lt 9 lw .5
END

