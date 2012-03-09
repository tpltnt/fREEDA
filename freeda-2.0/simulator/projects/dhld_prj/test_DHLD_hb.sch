<Qucs Schematic 0.1>
<Properties>
  <View=-20,-29,682,484,1,0,0>
  <Grid=10,10,1>
  <DataSet=test_DHLD_hb.dat>
  <DataDisplay=test_DHLD_hb.dpl>
  <OpenDisplay=1>
</Properties>
<Symbol>
</Symbol>
<Components>
  <vsource vsource:vcurrs 1 260 160 -26 18 0 0 "" 0 "0" 0 "0" 1 "0" 0 "0" 0 "0" 0 "0" 0>
  <r r:r1 1 390 230 15 -26 0 1 "" 0 "1e6" 1>
  <dhld dhld:d1 1 490 230 -26 34 0 0 "" 0 "2" 0 "0.468" 0 "2.54e-25" 0 "0.01813" 0 "6.92" 0 "2.25e-09" 0 "1e-11" 1 "1.65" 0 "1.79e-29" 0 "0.125" 0 "29.4" 0 "1.02e-13" 0 "1e+21" 0 "1e-3" 1>
  <r r:r2 1 570 210 -26 15 0 0 "" 0 "1e6" 1>
  <isource isource:i1 1 150 220 20 -26 0 1 "" 0 "150m" 1 "4m" 1 "1e9" 1 "-90" 1 "0" 0 "0" 0>
  <.out plot1 1 320 10 -35 14 0 0 "Mode=plot" 0 "term" 1 "Three" 1 "" 0 "" 0 "vf" 1 "" 0 "" 0 "" 0 "" 0 "" 0 "invfft 8 repeat" 1 "v_hb.freeda" 1>
  <svhb svhb1 1 150 10 0 33 0 0 "40" 1 "1e9" 1 "1" 0 "0" 0 "0" 0 "1" 0 "0" 0 "0" 0 "0" 0 "0" 0 "0" 0 "0" 0 "0" 0>
  <.out plot2 1 480 10 -35 14 0 0 "Mode=plot" 0 "term" 1 "Four" 1 "" 0 "" 0 "vf" 1 "" 0 "" 0 "" 0 "" 0 "" 0 "invfft 8 repeat" 1 "sn_hb.freeda" 1>
  <LRT * 1 570 260 0 0 0 0>
  <GND * 1 290 300 0 0 0 0>
</Components>
<Wires>
  <460 250 460 290 "" 0 0 0 "">
  <520 210 540 210 "Four" 560 180 8 "">
  <290 160 390 160 "" 0 0 0 "">
  <460 160 460 210 "" 0 0 0 "">
  <390 290 460 290 "" 0 0 0 "">
  <390 260 390 290 "" 0 0 0 "">
  <390 160 460 160 "Three" 460 130 44 "">
  <390 160 390 200 "" 0 0 0 "">
  <150 160 150 190 "" 0 0 0 "">
  <150 160 230 160 "" 0 0 0 "">
  <150 250 150 290 "" 0 0 0 "">
  <150 290 290 290 "" 0 0 0 "">
  <610 210 610 260 "oref" 640 240 40 "">
  <600 210 610 210 "" 0 0 0 "">
  <540 260 570 260 "" 0 0 0 "">
  <540 250 540 260 "" 0 0 0 "">
  <520 250 540 250 "" 0 0 0 "">
  <570 260 610 260 "" 0 0 0 "">
  <290 290 390 290 "" 0 0 0 "">
  <290 290 290 300 "" 0 0 0 "">
</Wires>
<Diagrams>
</Diagrams>
<Paintings>
</Paintings>
