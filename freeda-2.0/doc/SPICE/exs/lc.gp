#set terminal postscript
set terminal postscript portrait monochrome "Courier" 24
#set ylabel "Vout/Vin"
set ylabel "|Zout| (Ohm)"
set xlabel "Frequency (Hz)"
set nokey
set logscale x
plot "lcz.output" using "%f%f"  with linespoints
