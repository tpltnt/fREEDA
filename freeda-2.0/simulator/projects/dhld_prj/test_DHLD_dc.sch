<Qucs Schematic 0.1>
<Properties>
  <View=-20,-29,773,484,1,0,0>
  <Grid=10,10,1>
  <DataSet=test_DHLD_dc.dat>
  <DataDisplay=test_DHLD_dc.dpl>
  <OpenDisplay=1>
</Properties>
<Symbol>
</Symbol>
<Components>
  <vsource vsource:vs 1 60 200 18 -26 0 1 "" 0 "0m" 1 "0" 0 "0" 0 "0" 0 "0" 0 "0" 0>
  <vsource vsource:vp 1 200 150 -26 18 0 0 "" 0 "0" 0 "0" 1 "0" 0 "0" 0 "0" 0 "0" 0>
  <vccs vccs:v1 1 320 180 -26 34 0 0 "" 0 "-1" 1 "1e16" 1 "1e16" 1>
  <vsource vsource:vcurrs 1 400 150 -26 18 0 0 "" 0 "0" 0 "0" 1 "0" 0 "0" 0 "0" 0 "0" 0>
  <dhld dhld:d1 1 490 230 -26 34 0 0 "" 0 "2" 0 "0.468" 0 "2.54e-25" 0 "0.01813" 0 "6.92" 0 "2.25e-09" 0 "1e-11" 0 "1.65" 0 "1.79e-29" 0 "0.125" 0 "29.4" 0 "1.02e-13" 0 "1e+21" 0 "0.001" 0>
  <open open:o1 1 570 210 -20 0 0 0>
  <GND * 1 260 350 0 0 0 0>
  <dc dc1 1 50 30 0 33 0 0 "vs" 1 "91m" 1 "200m" 1 ".02m" 1 "0.0001" 0 "0" 0>
  <.out plot1 1 210 30 -35 14 0 0 "Mode=plot" 0 "term" 1 "Three" 1 "" 0 "" 0 "vt" 1 "" 0 "" 0 "" 0 "" 0 "" 0 "" 0 "DHLD_v_beta_5e_2_dc_freeda" 1>
  <.out plot2 1 470 30 -35 14 0 0 "Mode=plot" 0 "term" 1 "Four" 1 "" 0 "" 0 "vt" 1 "" 0 "" 0 "" 0 "" 0 "" 0 "" 0 "DHLD_sn_beat_5e_2_dc_freeda" 1>
  <LRT * 1 560 300 0 0 0 0>
</Components>
<Wires>
  <60 150 170 150 "" 0 0 0 "">
  <60 150 60 170 "" 0 0 0 "">
  <230 150 290 150 "" 0 0 0 "">
  <290 210 290 290 "" 0 0 0 "">
  <290 290 350 290 "" 0 0 0 "">
  <350 210 350 290 "" 0 0 0 "">
  <60 230 60 290 "" 0 0 0 "">
  <60 290 260 290 "" 0 0 0 "">
  <350 150 370 150 "" 0 0 0 "">
  <520 300 560 300 "" 0 0 0 "">
  <460 150 460 210 "Three" 490 160 40 "">
  <430 150 460 150 "" 0 0 0 "">
  <520 250 520 300 "" 0 0 0 "">
  <350 290 460 290 "" 0 0 0 "">
  <460 250 460 290 "" 0 0 0 "">
  <520 210 540 210 "Four" 560 180 8 "">
  <600 210 600 300 "oref" 640 240 64 "">
  <260 290 290 290 "" 0 0 0 "">
  <260 290 260 350 "" 0 0 0 "">
  <560 300 600 300 "" 0 0 0 "">
</Wires>
<Diagrams>
</Diagrams>
<Paintings>
</Paintings>
