reset
set terminal postscript landscape enhanced color "Times_Roman" 20
set encoding iso_8859_1
set pointsize 1.2
set size 0.95,0.95

set output "vdis.eps"

qunit=1.0
eunit=1.000000
funit=1.000000
p0 = 0.000000
pp0 = 0.241581
p1 = 0.241581
pp1 = 0.241581
p2 = 0.483161
pp2 = 0.139477
p3 = 0.622638

set xtics ( '0' qunit*0.000000, 'H' qunit*0.241581, '{/Symbol G}' qunit*0.483161, 'N' qunit*0.622638)

set key left top
set ylabel "Frequency (THz)"

plot [x=0:qunit*p3*1.0001] [funit*0.000000:funit*11.000000] \
'vline.dat' using (qunit*$1):(funit*$2) notitle w l lt 4, \
 'vdis.out' index 0 using (qunit*p0+qunit*$1):(funit*$5) notitle w l lt -1, \
 '' index 0 using (qunit*p0+qunit*$1):(funit*$6) notitle w l lt -1, \
 '' index 0 using (qunit*p0+qunit*$1):(funit*$7) notitle w l lt -1, \
 '' index 0 using (qunit*p0+qunit*$1):(funit*$8) notitle w l lt -1, \
 '' index 0 using (qunit*p0+qunit*$1):(funit*$9) notitle w l lt -1, \
 '' index 0 using (qunit*p0+qunit*$1):(funit*$10) notitle w l lt -1, \
 '' index 0 using (qunit*p0+qunit*$1):(funit*$11) notitle w l lt -1, \
 '' index 0 using (qunit*p0+qunit*$1):(funit*$12) notitle w l lt -1, \
 '' index 0 using (qunit*p0+qunit*$1):(funit*$13) notitle w l lt -1, \
 '' index 0 using (qunit*p0+qunit*$1):(funit*$14) notitle w l lt -1, \
 '' index 0 using (qunit*p0+qunit*$1):(funit*$15) notitle w l lt -1, \
 '' index 0 using (qunit*p0+qunit*$1):(funit*$16) notitle w l lt -1, \
 '' index 1 using (qunit*p1+qunit*$1):(funit*$5) notitle w l lt -1, \
 '' index 1 using (qunit*p1+qunit*$1):(funit*$6) notitle w l lt -1, \
 '' index 1 using (qunit*p1+qunit*$1):(funit*$7) notitle w l lt -1, \
 '' index 1 using (qunit*p1+qunit*$1):(funit*$8) notitle w l lt -1, \
 '' index 1 using (qunit*p1+qunit*$1):(funit*$9) notitle w l lt -1, \
 '' index 1 using (qunit*p1+qunit*$1):(funit*$10) notitle w l lt -1, \
 '' index 1 using (qunit*p1+qunit*$1):(funit*$11) notitle w l lt -1, \
 '' index 1 using (qunit*p1+qunit*$1):(funit*$12) notitle w l lt -1, \
 '' index 1 using (qunit*p1+qunit*$1):(funit*$13) notitle w l lt -1, \
 '' index 1 using (qunit*p1+qunit*$1):(funit*$14) notitle w l lt -1, \
 '' index 1 using (qunit*p1+qunit*$1):(funit*$15) notitle w l lt -1, \
 '' index 1 using (qunit*p1+qunit*$1):(funit*$16) notitle w l lt -1, \
 '' index 2 using (qunit*p2+qunit*$1):(funit*$5) notitle w l lt -1, \
 '' index 2 using (qunit*p2+qunit*$1):(funit*$6) notitle w l lt -1, \
 '' index 2 using (qunit*p2+qunit*$1):(funit*$7) notitle w l lt -1, \
 '' index 2 using (qunit*p2+qunit*$1):(funit*$8) notitle w l lt -1, \
 '' index 2 using (qunit*p2+qunit*$1):(funit*$9) notitle w l lt -1, \
 '' index 2 using (qunit*p2+qunit*$1):(funit*$10) notitle w l lt -1, \
 '' index 2 using (qunit*p2+qunit*$1):(funit*$11) notitle w l lt -1, \
 '' index 2 using (qunit*p2+qunit*$1):(funit*$12) notitle w l lt -1, \
 '' index 2 using (qunit*p2+qunit*$1):(funit*$13) notitle w l lt -1, \
 '' index 2 using (qunit*p2+qunit*$1):(funit*$14) notitle w l lt -1, \
 '' index 2 using (qunit*p2+qunit*$1):(funit*$15) notitle w l lt -1, \
 '' index 2 using (qunit*p2+qunit*$1):(funit*$16) notitle w l lt -1

#
#
