set term postscript color portrait enhanced 11
set output 'processBayes.ps'
set multiplot
set colors podo
set lmargin 4

set size .35,.25
set yrange [0:]
set style fill solid

bin1(x)=width*(floor(x/width)+.0)
norm1=1./155.   # normalizace na frekvence (1/pocetdat)

set origin .0,.7
set xtics 5.
width=1.
set boxwidth width*.9
set title 'Variance reduction (%)'
plot 'processBayes.dat' using (bin1($2*100)):(norm1/width) notitle smooth freq with boxes lc 7
set xtics autofreq
set autoscale x

set origin .35,.7
set xtics .5
width=.2
set boxwidth width*.9
set title 'M_0 (x10^{15} Nm)'
M0aprior1=1.7e15/1.e15
Mwsigma=0.1
a=(2./3./log(10.)/Mwsigma)**2
set xrange [0:5]
plot 'processBayes.dat' using (bin1($11/1.e15)):(norm1/width) notitle smooth freq with boxes lc 7,\
1./(M0aprior1*exp(.5/a)*sqrt(2.*pi/a))*exp(-0.5*(2./3.*log(x/M0aprior1)/log(10.)/Mwsigma)**2) notitle w l lc 3 lw 3
set xtics autofreq

set origin .7,.7
set xtics 2.
width=.5
set boxwidth width*.9
set title 'Mean stress drop (MPa)'
set xrange [1:3]
plot 'processBayes.dat' using (bin1($3)):(norm1/width) notitle smooth freq with boxes lc 7
set xtics autofreq
set autoscale x

set origin .0,.45
set xtics .05
width=.02
set boxwidth width*.9
set title 'Nucleation zone size (km^2)'
set xrange [0:.2]
plot 'processBayes.dat' using (bin1($5)):(norm1/width) notitle smooth freq with boxes lc 7,\
(x<.04*pi?1./.04/pi:0) notitle w l lc 3 lw 3
set xtics autofreq
set autoscale x

set origin 0.35,.45
set xtics .3
width=.05
set boxwidth width*.9
set title 'Mean nucleation overstress (MPa)'
set xrange [0:1.2]
plot 'processBayes.dat' using (bin1($10)):(norm1/width) notitle smooth freq with boxes lc 7,\
(x<1.?1.:0) notitle w l lc 3 lw 3
set xtics autofreq
set autoscale x

set origin .7,.45
set title 'Fracture energy density (kJ/m^2)'
set xtics 1.
width=.5
set boxwidth width*.9
set xrange [0:15]
plot 'processBayes.dat' using (bin1($3/$17/1.e3)):(norm1/width) notitle smooth freq with boxes lc 7
set xtics autofreq
set autoscale x

set origin .0,.2
set xtics .2
width=.05
set boxwidth width*.9
#set title 'Radiated energy density (kJ/m^2)'
set title 'Radiated energy (10^{9} J)'
set xrange [:1.]
plot 'processBayes.dat' using (bin1($7/1.e9)):(norm1/width) notitle smooth freq with boxes lc 7
set xtics autofreq
set autoscale x

set origin .35,.2
set xtics .1
width=.025
set boxwidth width*.9
set title 'Rad.energy-to-moment ratio (x10^{-5})'
set xrange [:.5]
plot 'processBayes.dat' using (bin1($7/$11*1e5)):(norm1/width) notitle smooth freq with boxes lc 7
set xtics autofreq
set autoscale x

set origin .7,.2
set title 'Radiation efficiency'
set xtics 0.02
width=.005
set boxwidth width*.9
plot 'processBayes.dat' using (bin1($8)):(norm1/width) notitle smooth freq with boxes lc 7
set xtics autofreq
set autoscale x

set origin .0,-.05
set xtics 10.
width=10.
set boxwidth width*.9
set title 'Rupture radius (m)'
set xrange [50:200.]
#set autoscale x
#set xtics auto
plot 'processBayes.dat' using (bin1(sqrt($17/3.1415926))):(norm1/width) notitle smooth freq with boxes lc 7
set xtics autofreq
set autoscale x

set origin .35,-.05
set xtics .2
width=.1
set boxwidth width*.9
set title 'Rupture duration (s)'
set xrange [0:2.]
#set autoscale x
#set xtics auto
plot 'processBayes.dat' using (bin1($4)):(norm1/width) notitle smooth freq with boxes lc 7
set xtics autofreq
set autoscale x
