<Qucs Schematic 0.1>
<Properties>
  <View=0,-49,800,800,1,0,49>
  <Grid=10,10,1>
  <DataSet=thermalMMIC1.dat>
  <DataDisplay=thermalMMIC1.dpl>
  <OpenDisplay=1>
</Properties>
<Symbol>
</Symbol>
<Components>
  <thermalheatsinkmmic1 thermalheatsinkmmic1:t1 1 240 170 -26 21 0 0 "" 0 "nsteps" 1 "deltat" 1 "temp" 1 "1" 1 "0.0004" 0 "0.0004" 0 "0.0004" 0 "0.00022" 0 "0.00018" 0 "0.00022" 0 "0.00018" 0 "46" 0 "5320" 0 "350" 0 "1" 0 "1.22" 0 "0" 0 "1" 0>
  <GND * 1 260 280 0 0 0 0>
  <vsource vsource:v1 1 420 180 18 -26 0 1 "" 0 "temp" 1 "0" 0 "0" 0 "0" 0 "0" 0 "0" 0>
  <isource isource:i1 1 40 190 20 -26 0 1 "" 0 "0.4" 1 "0" 0 "0" 0 "0" 0 "0" 0 "0" 0>
  <r r:r1 1 240 90 -26 15 0 0 "" 0 "1e5" 1>
  <.top top1 1 50 30 -10 14 0 0 ".options deltat=1us nsteps=101 temp=300." 1>
  <.out plot1 1 80 300 -35 14 0 0 "Mode=plot" 0 "element" 1 "" 0 "thermalheatsinkmmic1:t1" 1 "0" 1 "it" 1 "" 0 "" 0 "" 0 "" 0 "" 0 "" 0 "input.power" 1>
  <.out plot2 1 350 300 -35 14 0 0 "Mode=plot" 0 "element" 1 "" 0 "thermalheatsinkmmic1:t1" 1 "0" 1 "ut" 1 "" 0 "" 0 "" 0 "" 0 "" 0 "" 0 "output.temperaturec" 1>
  <svtr svtr1 1 450 30 0 33 0 0 "1" 0 "100us" 1 "deltat" 1 "8" 1 "200" 0 "1e-08" 0 "0" 0 "0" 0 "0" 0 "1" 0 "1" 0 "100" 0 "0" 0 "0" 0 "0" 0>
</Components>
<Wires>
  <270 90 340 90 "" 0 0 0 "">
  <420 90 420 150 "" 0 0 0 "">
  <420 210 420 280 "" 0 0 0 "">
  <260 280 420 280 "" 0 0 0 "">
  <40 90 40 160 "" 0 0 0 "">
  <40 90 140 90 "" 0 0 0 "">
  <40 220 40 280 "" 0 0 0 "">
  <40 280 260 280 "" 0 0 0 "">
  <140 170 210 170 "" 0 0 0 "">
  <140 90 210 90 "" 0 0 0 "">
  <140 90 140 170 "" 0 0 0 "">
  <340 90 420 90 "" 0 0 0 "">
  <340 90 340 170 "" 0 0 0 "">
  <270 170 340 170 "" 0 0 0 "">
</Wires>
<Diagrams>
</Diagrams>
<Paintings>
</Paintings>
