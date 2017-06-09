#################################################################################
# Gnuplot
#
#################################################################################
set terminal postscript eps enhanced color size 6.6,3.6
set output "energy-heating.eps"
#
# Set multiple plot layout
#============================================
set multiplot layout 1,1 rowsfirst
#
# Data file
#============================================
set datafile separator ","
#============================================
# Plot 1
#============================================
# Title 
set title "Energy - Heat"
# Grid settings 
#set grid xtics lc rgb "#CCCCCC" lw 0.2 lt 1 
#set grid ytics lc rgb "#CCCCCC" lw 0.2 lt 1 
# Key settings 
set key box opaque 
set border back 
# Axis 
set xlabel "Maximum Variation (Bohr)"
set ylabel "Energy (Ry)"
#set xtics 10
set ytics nomirror tc lt 1
set y2tics nomirror tc lt 1
# Circles 
# Plot 
plot \
'energy-heating.csv' using (1 * $1):(1 * $2)  title 'EoS'  with points axes x1y1