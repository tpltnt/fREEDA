<Qucs Schematic 0.1>
<Properties>
  <View=0,0,800,800,1,0,0>
  <Grid=10,10,1>
  <DataSet=ac.dat>
  <DataDisplay=ac.dpl>
  <OpenDisplay=1>
</Properties>
<Symbol>
</Symbol>
<Components>
  <GND * 1 260 310 0 0 0 0>
  <ac ac1 1 20 20 0 33 0 0 "0.01e9" 1 "20e9" 1 "100" 1 "1" 0>
  <vsource vsource:v1 1 80 190 18 -26 0 1 "" 0 "0" 0 "1" 1 "0" 0 "0" 0 "0" 0 "0" 0>
  <.out plot1 1 170 30 -35 14 0 0 "Mode=plot" 0 "term" 1 "2" 1 "" 0 "" 0 "vf" 1 "" 0 "" 0 "" 0 "" 0 "" 0 "mag" 1 "v11.out" 1>
  <.out plot2 1 360 20 -35 14 0 0 "Mode=plot" 0 "element" 1 "" 0 "vsource:v1" 1 "0" 1 "if" 1 "" 0 "" 0 "" 0 "" 0 "" 0 "mag db" 1 "iv1.out" 1>
  <c c:c1 1 350 190 -64 -26 0 3 "" 0 "1e-12" 1 "1e-08" 0 "0" 0>
  <r r:r1 1 450 190 -47 -26 0 3 "" 0 "2e3" 1>
  <vccs vccs:v1 1 240 190 -26 34 0 0 "" 0 "1e-3" 1 "1e3" 1 "2e3" 1>
</Components>
<Wires>
  <80 160 210 160 "" 0 0 0 "">
  <270 160 350 160 "2" 310 130 13 "">
  <270 220 350 220 "" 0 0 0 "">
  <350 220 450 220 "" 0 0 0 "">
  <260 310 350 310 "" 0 0 0 "">
  <350 220 350 310 "" 0 0 0 "">
  <140 310 260 310 "" 0 0 0 "">
  <80 220 140 220 "" 0 0 0 "">
  <140 220 210 220 "" 0 0 0 "">
  <140 220 140 310 "" 0 0 0 "">
  <350 160 450 160 "" 0 0 0 "">
</Wires>
<Diagrams>
</Diagrams>
<Paintings>
</Paintings>