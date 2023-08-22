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
pp0 = 0.109277
p1 = 0.109277
pp1 = 0.128153
p2 = 0.237430
pp2 = 0.082400
p3 = 0.319830
pp3 = 0.061026
p4 = 0.380856

set xtics ( '{/Symbol G}' qunit*0.000000, 'X' qunit*0.109277, 'Y' qunit*0.237430, '{/Symbol G}' qunit*0.319830, 'Z' qunit*0.380856)

set key left top
set ylabel "Frequency (THz)"

plot [x=0:qunit*p4*1.0001] [funit*0.000000:funit*13.000000] \
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
 '' index 0 using (qunit*p0+qunit*$1):(funit*$17) notitle w l lt -1, \
 '' index 0 using (qunit*p0+qunit*$1):(funit*$18) notitle w l lt -1, \
 '' index 0 using (qunit*p0+qunit*$1):(funit*$19) notitle w l lt -1, \
 '' index 0 using (qunit*p0+qunit*$1):(funit*$20) notitle w l lt -1, \
 '' index 0 using (qunit*p0+qunit*$1):(funit*$21) notitle w l lt -1, \
 '' index 0 using (qunit*p0+qunit*$1):(funit*$22) notitle w l lt -1, \
 '' index 0 using (qunit*p0+qunit*$1):(funit*$23) notitle w l lt -1, \
 '' index 0 using (qunit*p0+qunit*$1):(funit*$24) notitle w l lt -1, \
 '' index 0 using (qunit*p0+qunit*$1):(funit*$25) notitle w l lt -1, \
 '' index 0 using (qunit*p0+qunit*$1):(funit*$26) notitle w l lt -1, \
 '' index 0 using (qunit*p0+qunit*$1):(funit*$27) notitle w l lt -1, \
 '' index 0 using (qunit*p0+qunit*$1):(funit*$28) notitle w l lt -1, \
 '' index 0 using (qunit*p0+qunit*$1):(funit*$29) notitle w l lt -1, \
 '' index 0 using (qunit*p0+qunit*$1):(funit*$30) notitle w l lt -1, \
 '' index 0 using (qunit*p0+qunit*$1):(funit*$31) notitle w l lt -1, \
 '' index 0 using (qunit*p0+qunit*$1):(funit*$32) notitle w l lt -1, \
 '' index 0 using (qunit*p0+qunit*$1):(funit*$33) notitle w l lt -1, \
 '' index 0 using (qunit*p0+qunit*$1):(funit*$34) notitle w l lt -1, \
 '' index 0 using (qunit*p0+qunit*$1):(funit*$35) notitle w l lt -1, \
 '' index 0 using (qunit*p0+qunit*$1):(funit*$36) notitle w l lt -1, \
 '' index 0 using (qunit*p0+qunit*$1):(funit*$37) notitle w l lt -1, \
 '' index 0 using (qunit*p0+qunit*$1):(funit*$38) notitle w l lt -1, \
 '' index 0 using (qunit*p0+qunit*$1):(funit*$39) notitle w l lt -1, \
 '' index 0 using (qunit*p0+qunit*$1):(funit*$40) notitle w l lt -1, \
 '' index 0 using (qunit*p0+qunit*$1):(funit*$41) notitle w l lt -1, \
 '' index 0 using (qunit*p0+qunit*$1):(funit*$42) notitle w l lt -1, \
 '' index 0 using (qunit*p0+qunit*$1):(funit*$43) notitle w l lt -1, \
 '' index 0 using (qunit*p0+qunit*$1):(funit*$44) notitle w l lt -1, \
 '' index 0 using (qunit*p0+qunit*$1):(funit*$45) notitle w l lt -1, \
 '' index 0 using (qunit*p0+qunit*$1):(funit*$46) notitle w l lt -1, \
 '' index 0 using (qunit*p0+qunit*$1):(funit*$47) notitle w l lt -1, \
 '' index 0 using (qunit*p0+qunit*$1):(funit*$48) notitle w l lt -1, \
 '' index 0 using (qunit*p0+qunit*$1):(funit*$49) notitle w l lt -1, \
 '' index 0 using (qunit*p0+qunit*$1):(funit*$50) notitle w l lt -1, \
 '' index 0 using (qunit*p0+qunit*$1):(funit*$51) notitle w l lt -1, \
 '' index 0 using (qunit*p0+qunit*$1):(funit*$52) notitle w l lt -1, \
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
 '' index 1 using (qunit*p1+qunit*$1):(funit*$17) notitle w l lt -1, \
 '' index 1 using (qunit*p1+qunit*$1):(funit*$18) notitle w l lt -1, \
 '' index 1 using (qunit*p1+qunit*$1):(funit*$19) notitle w l lt -1, \
 '' index 1 using (qunit*p1+qunit*$1):(funit*$20) notitle w l lt -1, \
 '' index 1 using (qunit*p1+qunit*$1):(funit*$21) notitle w l lt -1, \
 '' index 1 using (qunit*p1+qunit*$1):(funit*$22) notitle w l lt -1, \
 '' index 1 using (qunit*p1+qunit*$1):(funit*$23) notitle w l lt -1, \
 '' index 1 using (qunit*p1+qunit*$1):(funit*$24) notitle w l lt -1, \
 '' index 1 using (qunit*p1+qunit*$1):(funit*$25) notitle w l lt -1, \
 '' index 1 using (qunit*p1+qunit*$1):(funit*$26) notitle w l lt -1, \
 '' index 1 using (qunit*p1+qunit*$1):(funit*$27) notitle w l lt -1, \
 '' index 1 using (qunit*p1+qunit*$1):(funit*$28) notitle w l lt -1, \
 '' index 1 using (qunit*p1+qunit*$1):(funit*$29) notitle w l lt -1, \
 '' index 1 using (qunit*p1+qunit*$1):(funit*$30) notitle w l lt -1, \
 '' index 1 using (qunit*p1+qunit*$1):(funit*$31) notitle w l lt -1, \
 '' index 1 using (qunit*p1+qunit*$1):(funit*$32) notitle w l lt -1, \
 '' index 1 using (qunit*p1+qunit*$1):(funit*$33) notitle w l lt -1, \
 '' index 1 using (qunit*p1+qunit*$1):(funit*$34) notitle w l lt -1, \
 '' index 1 using (qunit*p1+qunit*$1):(funit*$35) notitle w l lt -1, \
 '' index 1 using (qunit*p1+qunit*$1):(funit*$36) notitle w l lt -1, \
 '' index 1 using (qunit*p1+qunit*$1):(funit*$37) notitle w l lt -1, \
 '' index 1 using (qunit*p1+qunit*$1):(funit*$38) notitle w l lt -1, \
 '' index 1 using (qunit*p1+qunit*$1):(funit*$39) notitle w l lt -1, \
 '' index 1 using (qunit*p1+qunit*$1):(funit*$40) notitle w l lt -1, \
 '' index 1 using (qunit*p1+qunit*$1):(funit*$41) notitle w l lt -1, \
 '' index 1 using (qunit*p1+qunit*$1):(funit*$42) notitle w l lt -1, \
 '' index 1 using (qunit*p1+qunit*$1):(funit*$43) notitle w l lt -1, \
 '' index 1 using (qunit*p1+qunit*$1):(funit*$44) notitle w l lt -1, \
 '' index 1 using (qunit*p1+qunit*$1):(funit*$45) notitle w l lt -1, \
 '' index 1 using (qunit*p1+qunit*$1):(funit*$46) notitle w l lt -1, \
 '' index 1 using (qunit*p1+qunit*$1):(funit*$47) notitle w l lt -1, \
 '' index 1 using (qunit*p1+qunit*$1):(funit*$48) notitle w l lt -1, \
 '' index 1 using (qunit*p1+qunit*$1):(funit*$49) notitle w l lt -1, \
 '' index 1 using (qunit*p1+qunit*$1):(funit*$50) notitle w l lt -1, \
 '' index 1 using (qunit*p1+qunit*$1):(funit*$51) notitle w l lt -1, \
 '' index 1 using (qunit*p1+qunit*$1):(funit*$52) notitle w l lt -1, \
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
 '' index 2 using (qunit*p2+qunit*$1):(funit*$16) notitle w l lt -1, \
 '' index 2 using (qunit*p2+qunit*$1):(funit*$17) notitle w l lt -1, \
 '' index 2 using (qunit*p2+qunit*$1):(funit*$18) notitle w l lt -1, \
 '' index 2 using (qunit*p2+qunit*$1):(funit*$19) notitle w l lt -1, \
 '' index 2 using (qunit*p2+qunit*$1):(funit*$20) notitle w l lt -1, \
 '' index 2 using (qunit*p2+qunit*$1):(funit*$21) notitle w l lt -1, \
 '' index 2 using (qunit*p2+qunit*$1):(funit*$22) notitle w l lt -1, \
 '' index 2 using (qunit*p2+qunit*$1):(funit*$23) notitle w l lt -1, \
 '' index 2 using (qunit*p2+qunit*$1):(funit*$24) notitle w l lt -1, \
 '' index 2 using (qunit*p2+qunit*$1):(funit*$25) notitle w l lt -1, \
 '' index 2 using (qunit*p2+qunit*$1):(funit*$26) notitle w l lt -1, \
 '' index 2 using (qunit*p2+qunit*$1):(funit*$27) notitle w l lt -1, \
 '' index 2 using (qunit*p2+qunit*$1):(funit*$28) notitle w l lt -1, \
 '' index 2 using (qunit*p2+qunit*$1):(funit*$29) notitle w l lt -1, \
 '' index 2 using (qunit*p2+qunit*$1):(funit*$30) notitle w l lt -1, \
 '' index 2 using (qunit*p2+qunit*$1):(funit*$31) notitle w l lt -1, \
 '' index 2 using (qunit*p2+qunit*$1):(funit*$32) notitle w l lt -1, \
 '' index 2 using (qunit*p2+qunit*$1):(funit*$33) notitle w l lt -1, \
 '' index 2 using (qunit*p2+qunit*$1):(funit*$34) notitle w l lt -1, \
 '' index 2 using (qunit*p2+qunit*$1):(funit*$35) notitle w l lt -1, \
 '' index 2 using (qunit*p2+qunit*$1):(funit*$36) notitle w l lt -1, \
 '' index 2 using (qunit*p2+qunit*$1):(funit*$37) notitle w l lt -1, \
 '' index 2 using (qunit*p2+qunit*$1):(funit*$38) notitle w l lt -1, \
 '' index 2 using (qunit*p2+qunit*$1):(funit*$39) notitle w l lt -1, \
 '' index 2 using (qunit*p2+qunit*$1):(funit*$40) notitle w l lt -1, \
 '' index 2 using (qunit*p2+qunit*$1):(funit*$41) notitle w l lt -1, \
 '' index 2 using (qunit*p2+qunit*$1):(funit*$42) notitle w l lt -1, \
 '' index 2 using (qunit*p2+qunit*$1):(funit*$43) notitle w l lt -1, \
 '' index 2 using (qunit*p2+qunit*$1):(funit*$44) notitle w l lt -1, \
 '' index 2 using (qunit*p2+qunit*$1):(funit*$45) notitle w l lt -1, \
 '' index 2 using (qunit*p2+qunit*$1):(funit*$46) notitle w l lt -1, \
 '' index 2 using (qunit*p2+qunit*$1):(funit*$47) notitle w l lt -1, \
 '' index 2 using (qunit*p2+qunit*$1):(funit*$48) notitle w l lt -1, \
 '' index 2 using (qunit*p2+qunit*$1):(funit*$49) notitle w l lt -1, \
 '' index 2 using (qunit*p2+qunit*$1):(funit*$50) notitle w l lt -1, \
 '' index 2 using (qunit*p2+qunit*$1):(funit*$51) notitle w l lt -1, \
 '' index 2 using (qunit*p2+qunit*$1):(funit*$52) notitle w l lt -1, \
 '' index 3 using (qunit*p3+qunit*$1):(funit*$5) notitle w l lt -1, \
 '' index 3 using (qunit*p3+qunit*$1):(funit*$6) notitle w l lt -1, \
 '' index 3 using (qunit*p3+qunit*$1):(funit*$7) notitle w l lt -1, \
 '' index 3 using (qunit*p3+qunit*$1):(funit*$8) notitle w l lt -1, \
 '' index 3 using (qunit*p3+qunit*$1):(funit*$9) notitle w l lt -1, \
 '' index 3 using (qunit*p3+qunit*$1):(funit*$10) notitle w l lt -1, \
 '' index 3 using (qunit*p3+qunit*$1):(funit*$11) notitle w l lt -1, \
 '' index 3 using (qunit*p3+qunit*$1):(funit*$12) notitle w l lt -1, \
 '' index 3 using (qunit*p3+qunit*$1):(funit*$13) notitle w l lt -1, \
 '' index 3 using (qunit*p3+qunit*$1):(funit*$14) notitle w l lt -1, \
 '' index 3 using (qunit*p3+qunit*$1):(funit*$15) notitle w l lt -1, \
 '' index 3 using (qunit*p3+qunit*$1):(funit*$16) notitle w l lt -1, \
 '' index 3 using (qunit*p3+qunit*$1):(funit*$17) notitle w l lt -1, \
 '' index 3 using (qunit*p3+qunit*$1):(funit*$18) notitle w l lt -1, \
 '' index 3 using (qunit*p3+qunit*$1):(funit*$19) notitle w l lt -1, \
 '' index 3 using (qunit*p3+qunit*$1):(funit*$20) notitle w l lt -1, \
 '' index 3 using (qunit*p3+qunit*$1):(funit*$21) notitle w l lt -1, \
 '' index 3 using (qunit*p3+qunit*$1):(funit*$22) notitle w l lt -1, \
 '' index 3 using (qunit*p3+qunit*$1):(funit*$23) notitle w l lt -1, \
 '' index 3 using (qunit*p3+qunit*$1):(funit*$24) notitle w l lt -1, \
 '' index 3 using (qunit*p3+qunit*$1):(funit*$25) notitle w l lt -1, \
 '' index 3 using (qunit*p3+qunit*$1):(funit*$26) notitle w l lt -1, \
 '' index 3 using (qunit*p3+qunit*$1):(funit*$27) notitle w l lt -1, \
 '' index 3 using (qunit*p3+qunit*$1):(funit*$28) notitle w l lt -1, \
 '' index 3 using (qunit*p3+qunit*$1):(funit*$29) notitle w l lt -1, \
 '' index 3 using (qunit*p3+qunit*$1):(funit*$30) notitle w l lt -1, \
 '' index 3 using (qunit*p3+qunit*$1):(funit*$31) notitle w l lt -1, \
 '' index 3 using (qunit*p3+qunit*$1):(funit*$32) notitle w l lt -1, \
 '' index 3 using (qunit*p3+qunit*$1):(funit*$33) notitle w l lt -1, \
 '' index 3 using (qunit*p3+qunit*$1):(funit*$34) notitle w l lt -1, \
 '' index 3 using (qunit*p3+qunit*$1):(funit*$35) notitle w l lt -1, \
 '' index 3 using (qunit*p3+qunit*$1):(funit*$36) notitle w l lt -1, \
 '' index 3 using (qunit*p3+qunit*$1):(funit*$37) notitle w l lt -1, \
 '' index 3 using (qunit*p3+qunit*$1):(funit*$38) notitle w l lt -1, \
 '' index 3 using (qunit*p3+qunit*$1):(funit*$39) notitle w l lt -1, \
 '' index 3 using (qunit*p3+qunit*$1):(funit*$40) notitle w l lt -1, \
 '' index 3 using (qunit*p3+qunit*$1):(funit*$41) notitle w l lt -1, \
 '' index 3 using (qunit*p3+qunit*$1):(funit*$42) notitle w l lt -1, \
 '' index 3 using (qunit*p3+qunit*$1):(funit*$43) notitle w l lt -1, \
 '' index 3 using (qunit*p3+qunit*$1):(funit*$44) notitle w l lt -1, \
 '' index 3 using (qunit*p3+qunit*$1):(funit*$45) notitle w l lt -1, \
 '' index 3 using (qunit*p3+qunit*$1):(funit*$46) notitle w l lt -1, \
 '' index 3 using (qunit*p3+qunit*$1):(funit*$47) notitle w l lt -1, \
 '' index 3 using (qunit*p3+qunit*$1):(funit*$48) notitle w l lt -1, \
 '' index 3 using (qunit*p3+qunit*$1):(funit*$49) notitle w l lt -1, \
 '' index 3 using (qunit*p3+qunit*$1):(funit*$50) notitle w l lt -1, \
 '' index 3 using (qunit*p3+qunit*$1):(funit*$51) notitle w l lt -1, \
 '' index 3 using (qunit*p3+qunit*$1):(funit*$52) notitle w l lt -1

#
#
