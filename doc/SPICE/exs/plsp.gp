#set terminal postscript
set terminal postscript portrait monochrome "Courier" 24
set xlabel "Time (ns)"
set ylabel "Voltage (V)"
set nokey
set label 1 "V(1)" at 11,1.4
set label 2 "V(2)" at 14,0.8
plot "simple.out2" using "%f%f"  with linespoints, \
"simple.out2" using "%f%*f%f"  with linespoints
