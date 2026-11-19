scan = int(ARG1)
dir = ARG2

set title sprintf("Convergence of absolute translation error in ICP algorithms\nScan %03d", scan)
set grid
set xlabel 'Iterations'
set ylabel 'Error [mm]'
set xrange [0:50]
set xtics 10
set terminal jpeg

set style line 1 \
    linecolor rgb '#0060ad' \
    linetype 1 linewidth 2 \
    pointtype 7 pointsize 1.5
set style line 2 \
    linecolor rgb '#dd181f' \
    linetype 1 linewidth 2 \
    pointtype 5 pointsize 1.5

set output sprintf("%s/t_scan%03d.jpeg", dir, scan)
    
plot sprintf("%s/fix_scan%03d_trans.dat", dir, scan)\
	using 1:2 with linespoints ls 1\
	title "Fixed-point", \
	sprintf("%s/ref_scan%03d_trans.dat", dir, scan)\
	using 1:2 with linespoints ls 2 \
	title "Floating-Point"
