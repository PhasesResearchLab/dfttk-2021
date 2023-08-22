reset
set terminal png font Times_Roman 96 size 4096,3072 linewidth 7
set encoding utf8
set pointsize 1.2
set xlabel "Frequency (THz)"
set ylabel "Phonon DOS (1/THz)"
set size 0.95,0.95
set key left top

set output "vdos.png"

set xlabel "Phonon frequncy (THz)"
set ylabel "Phonon DOS (1/THz/cell)"

funit=1.000000
plot 'vdos.out' using ($1*funit*1.e-12):($2/funit*1.e12) notitle w l lt -1

