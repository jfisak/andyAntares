set terminal postscript eps color dashed dashlength 5 25
set output "dynamicalCells.eps"
set encoding iso_8859_1
set size  1,1
set view equal xyz
set xtics 4e11
set ytics 4e11
set ztics 4e11

set style line 1 lc rgb 'green'
set style arrow 1 head filled ls 1 

splot 'photpos_dyn1.dat' us 7:8:9 notitle,\
	'photpos_dyn1.dat' us 1:2:3:4:(0):(0) with vectors arrowstyle 1 notitle,\
	'photpos_dyn1.dat' us 1:2:3:(0):5:(0) with vectors arrowstyle 1 notitle,\
	'photpos_dyn1.dat' us 1:2:3:(0):(0):6 with vectors arrowstyle 1 notitle
