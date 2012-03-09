#set terminal postscript
set terminal postscript portrait monochrome "Courier" 24
set xlabel "Vout (V)"
set ylabel "Output Current (A)"
set nokey
plot "cmos_inv.out" using "%f%f"  with lines
