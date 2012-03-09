lspiral readme
2008.04.16
George Feller
grfeller@ncsu.edu

lspiral.doc - documenatation about the lsprial model


A little bit about the other files contained in the document directory:

A. Netlists:

lsp_ac_y_11.net: Admittance parameter extraction (Y11) to determine R11 and L11 vs. frequency.

lsp_ac_y_12.net: Admittance parameter extraction (Y12) to determine R12 and L12 vs. frequency.

lsp_proto.net: Lumped circuit version of lspiral.

lsp.net: Transient response 

lsp_ideal.net: Ideal transient response

B. Data files:

lsp_tr.dat - transient response of RLC tank using lspira
lsp_tr_ideal.dat - transient response of ideal RLC tan
y11.dat - Y11 admittance data
y12.dat - Y12 admittance data

C. gnuplot command files: These gnuplot commands side-step some minor plotting quirks exhibited by freeda (most notably, the automatic conversion of all characters to lower case, which leads to unsatisfying labels).

usage: gnuplot <filename>

e.g. gnuplot plot_lt_spiral.cmd

plot_l11.cmd - inductance computed from admittance Y11
plot_l12.cmd - inductance computed from admittance Y12
plot_lsp_tr_ideal.cmd - transient response of ideal RLC tank
plot_lsp_tr.cmd - transient response of RLC tank using lspiral
plot_q11.cmd - quality factor (q11) from admittance Y11
plot_r11.cmd - resistance from admittance Y11
plot_r12.cmd - resistance from admittance Y12

D. Postscript plot files (as above)
l11.ps
l12.ps
lsp_lt.ps - same as plot_lt_spiral above
lt_ideal.ps - same as ideal RLC tank above
q11.ps
r11.ps
r12.ps


