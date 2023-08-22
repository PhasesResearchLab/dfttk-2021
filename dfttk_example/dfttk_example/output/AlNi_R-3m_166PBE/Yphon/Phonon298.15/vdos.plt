reset
set terminal postscript landscape enhanced color "Times_Roman" 20
set encoding iso_8859_1
set pointsize 1.2
set xlabel "Frequency (THz)"
set ylabel "Phonon DOS (1/THz)"
set size 0.95,0.95
set key left top

set output "vdos.eps"

funit=1.000000
plot 'vdos.out' using ($1*funit*1.e-12):($2/funit*1.e12) notitle w l lt -1

