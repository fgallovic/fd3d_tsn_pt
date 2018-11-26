#!/bin/bash
nxt=300
nzt=140  #po odecteni 2 kvuli povrchu
dh=0.1
L=30.
W=14.

python << END
from numpy import *
import matplotlib.pyplot as plt
slip=loadtxt('processBayes.slipmodels.dat')
stex=loadtxt('processBayes.strengthexcess.dat')
NM=slip.size/$nxt/$nzt
slip=slip.reshape(($nzt,$nxt,NM))
stex=stex.reshape(($nzt,$nxt,NM))
meanslip=slip.mean(axis=2)

x=linspace($dh/2.,$L-$dh/2.,$nxt)
z=linspace($dh/2.,$W-$dh/2.,$nzt)
f= open('processBayes.slipcontour.dat','w')
f2=open('processBayes.stexcontour.dat','w')
for i in range(NM):
  slip1=slip[:,:,i]
  stex1=stex[:,:,i]
  slipmax=slip1.max()
  cnt=plt.contour(x,z,slip1,[slipmax*.05]).collections[0].get_paths()[0].vertices
  cnt1=plt.contour(x,z,stex1,[0.]).collections[0].get_paths()[0].vertices
  f.write('\n'.join([str(s[0])+' '+str(s[1]) for s in cnt]))
  f.write('\n\n')
  f2.write('\n'.join([str(s[0])+' '+str(s[1]) for s in cnt1]))
  f2.write('\n\n')
f2.close()
f.write('\n')
f.write('\n')
slipmax=meanslip.max()
cnt=plt.contour(x,z,meanslip,[slipmax*.05]).collections[0].get_paths()[0].vertices
f.write('\n'.join([str(s[0])+' '+str(s[1]) for s in cnt]))
f.close()
END

gnuplot5 << END
set term postscript color solid enhanced 12
set output 'processBayes.slipmodels.ps'
#set term png truecolor enhanced 12
#set output 'processBayes.slipmodels.png'
set multiplot
set size 0.5,0.5
set palette defined ( 0 "white", 2 "skyblue", 3 "light-green", 6 "yellow", 10 "light-red" )
set xtics out scale .2
set ytics out scale .2
set size ratio -1
#set cbrange [0:$MAX]
set xrange [0:$nxt*$dh]
set yrange [0:$nzt*$dh]
set xlabel 'Along strike (km)'
set ylabel 'Up dip (km)'
set cbrange [0:]

set origin 0.,0.5
set title 'Mean slip (m)'
plot 'processBayes.meansigma2.dat' matrix u (\$1*$dh):(\$2*$dh):3 index 0 notitle w image,\
'processBayes.slipcontour.dat' index 0 w l notitle lt 9 lw .1 lc rgb "#EE000000",\
'processBayes.slipcontour.dat' index 1 w l notitle lt -1 lw 1
#'processBayes.stexcontour.dat' index 0 w l notitle lt 9 lw .1 lc 3,\  #nejak nefunguje

END
