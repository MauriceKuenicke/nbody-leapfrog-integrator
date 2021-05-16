#!/bin/bash
#
# shell-script to plot several snapshots from a single data file
# the data is in a file called 'out' and has the format
#
#     time    x1  y1  z1     x2  y2  z2     x3  y3  z3
#
# settings allow you to change the name of the input file and also
# skip lines for plotting (e.g. every=10 plots only one in ten lines)
# the resulting plots are in PNG-format and store in numbered files
# a movie can be created with the make_movie.sh script

# settings
infile='out.dat'
every=1000

# loop over lines
i=0     # file number counter

cat << EOF > tmp.awk
{
 if (NR%$every==0) print
}
EOF

awk -f tmp.awk $infile |\
while read line
do

    # set name of output png file
    outfile=`printf "plots/threebody_%5.5i.png" $i`
    echo generating $outfile


    echo $line > tmp.dat

# calling gnuplot to make plot
# use lots of settings in gnuplot including multiplot
# make sure FreeSans.ttf is in the stated dir or change accordingly
gnuplot << EOF
set term png font "/usr/share/fonts/truetype/freefont/FreeSans.ttf" 20 size 1024,480
set out '$outfile'

set multiplot

unset key
set border 31 lw 3
set size square

set origin 0.0 , 0.0
set size 0.5, 1

set xlabel 'x'
set ylabel 'y'

set xrange [-2 to 2]
set yrange [-2 to 2]

set xtics 1
set ytics 1


plot 'tmp.dat' using 2:3 with linespoints pt 7 ps 3,\
     'tmp.dat' using 5:6 with linespoints pt 7 ps 3,\
     'tmp.dat' using 8:9 with linespoints pt 7 ps 3

set origin 0.45 , 0.0
set size 0.6, 1

set xrange [-0.5 to 0.5]
set yrange [-0.5 to 0.5]

set xtics 0.5
set ytics 0.5

plot 'tmp.dat' using 2:3 with linespoints pt 7 ps 3,\
     'tmp.dat' using 5:6 with linespoints pt 7 ps 3,\
     'tmp.dat' using 8:9 with linespoints pt 7 ps 3

EOF

((i++))
done

rm tmp.dat tmp.awk