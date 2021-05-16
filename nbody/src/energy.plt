# Plot with gnuplot -p energy.plt
# Plots the relative energy error of the nbody system over the timesteps
# Only use every 10000th data point --> can be changed in line 15

set terminal png
set xlabel "Timestep"
set ylabel "Partial Energy Error"
set output "plots/energy_difference.png"
m="out.dat"

set nokey

set yrange[0 to 0.00000000000002]

plot m every 100000 using 1:12 with lines lw 2 title "Partial Energy Error"
