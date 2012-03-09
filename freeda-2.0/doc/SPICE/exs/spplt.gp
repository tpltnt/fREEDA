#set terminal postscript
set terminal postscript portrait monochrome "Courier" 24
set xlabel "Time (s)"
set ylabel "Cost"
set nokey
set title "Partition Cost vs. Number Containers"

set xrange [2:15]
set label 1 "partition cost, level 2" at 11.2,137
set label 2 "seeding estimate" at 11.2,220
set label 3 "increm. seeding estimate" at 11.2,145
set label 4 "partition cost, level 1" at 10.2,150
set output "basic4plot.ps"
set title "Partition Cost vs. Number Containers, test basic4"
plot "basic4_kl_lev2.dat" with linespoints, "basic4_kl_lev1.dat" with linespoints, "basic4_increm.dat" with linespoints, "basic4_seeding.dat" with linespoints

set label 1 at 11.2,115
set label 2 at 11.2,165
set label 3 at 11.2,95
set label 4 at 8.2,106
set output "basic5plot.ps"
set title "Partition Cost vs. Number Containers, test 5"
plot "basic5_kl_lev2.dat" with linespoints,"basic5_kl_lev1.dat" with linespoints,"basic5_increm.dat" with linespoints,"basic5_seeding.dat" with linespoints

set label 1 at 11.2,115
set label 2 at 11.2,155
set label 3 at 11.2,90
set label 4 at 8.2,106
set output "basic6plot.ps"
set title "Partition Cost vs. Number Containers, test 6"
plot "basic6_kl_lev2.dat" with linespoints,"basic6_kl_lev1.dat" with linespoints,"basic6_increm.dat" with linespoints,"basic6_seeding.dat" with linespoints


set label 1 at 11.2,110
set label 2 at 11.2,165
set label 3 at 11.2,100
set label 4 at 8.2,114
set output "basic7plot.ps"
set title "Partition Cost vs. Number Containers, test 7"
plot "basic7_kl_lev2.dat" with linespoints,"basic7_kl_lev1.dat" with linespoints,"basic7_increm.dat" with linespoints,"basic7_seeding.dat" with linespoints

