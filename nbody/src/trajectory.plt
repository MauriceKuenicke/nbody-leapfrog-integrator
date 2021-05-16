# Plot with gnuplot -p trajectory.plt
# Plots the tajectories of the 3 bodies in the system

set terminal png
set xlabel "x"
set ylabel "y"
set output "plots/trajectories.png"
m="out.dat"

set xrange[-2 to 2]
set yrange[-2 to 2]

set style data lines
set key top left

plot m  every 100 using 2:3 lw 2 lt rgb "orange" title "Body 1",\
     m  every 100 using 5:6 lw 2 lt rgb "blue" title "Body 2",\
     m  every 100 using 8:9 lw 2 lt rgb "green" title "Body 3",\
     "< tail -1 out.dat" using 2:3 with points lt rgb "orange" pt 7 ps 2 notitle,\
     "< tail -1 out.dat" using 5:6 with points lt rgb "blue" pt 7 ps 2 notitle,\
     "< tail -1 out.dat" using 8:9 with points lt rgb "green" pt 7 ps 2 notitle
