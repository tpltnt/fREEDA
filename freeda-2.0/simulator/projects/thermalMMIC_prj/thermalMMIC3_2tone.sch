<Qucs Schematic 0.1>
<Properties>
  <View=-4,-102,883,600,1,56,198>
  <Grid=10,10,1>
  <DataSet=thermalMMIC3_2tone.dat>
  <DataDisplay=thermalMMIC3_2tone.dpl>
  <OpenDisplay=1>
</Properties>
<Symbol>
</Symbol>
<Components>
  <r r:rin 1 180 200 -26 15 0 0 "" 0 "50" 1>
  <GND * 1 260 340 0 0 0 0>
  <vsource vsource:vbias 1 200 110 18 -26 0 1 "" 0 "-0.2" 1 "0" 0 "0" 0 "0" 0 "0" 0 "0" 0>
  <GND * 1 200 140 0 0 0 0>
  <r r:r2 1 260 60 -26 -45 0 2 "" 0 "100" 1>
  <r r:r3 1 470 10 -40 -26 0 3 "" 0 "10" 1>
  <c c:cload 1 540 100 -26 17 0 0 "" 0 "20e-12" 1 "1e-08" 0 "0" 0>
  <r r:rload 1 620 100 -26 15 0 0 "" 0 "50" 1>
  <GND * 1 670 100 0 0 0 0>
  <vsource vsource:vdrain 1 530 20 18 -26 0 1 "" 0 "3" 1 "0" 0 "0" 0 "0" 0 "0" 0 "0" 0>
  <GND * 1 530 50 0 0 0 0>
  <r r:rs 1 470 310 -57 -26 0 3 "" 0 "1.144" 1>
  <vsource vsource:t1 1 770 220 18 -26 0 1 "" 0 "temp" 1 "0" 0 "0" 0 "0" 0 "0" 0 "0" 0>
  <thermalheatsinkmmic1 thermalheatsinkmmic1:t1 1 630 180 -26 21 0 0 "" 0 "0" 0 "0" 0 "temp" 1 "0" 0 "0.0004" 0 "0.0004" 0 "0.0004" 0 "0.00022" 0 "0.00018" 0 "0.00022" 0 "0.00018" 0 "46" 0 "5320" 0 "350" 0 "50" 1 "1.22" 0 "0" 0 "1" 0>
  <mesfetct mesfetct:m1 1 470 230 30 -26 0 0 "mesfet" 1 "0.016542" 0 "0.0500214" 0 "0.02012" 0 "-0.00806" 0 "-0.0394707" 0 "4" 0 "2.16505" 0 "-1e+10" 0 "0.52785e-12" 0 "0.087e-12" 0 "1e-9" 0 "1" 0 "1e-9" 0 "10" 0 "5e-12" 0 "0.8" 0 "0.5" 0 "15" 0 "293" 0 "0" 0 "0" 0 "0" 0 "0" 0 "0" 0 "0.8" 0 "0.5" 0 "2" 0 "0" 0 "1.22" 0 "300" 0 "1" 0>
  <.model model1 1 530 290 -29 14 0 0 "mesfet.mdl" 1>
  <l l:l2 1 320 60 -26 -40 0 2 "" 0 "15e-9" 1 "1e-08" 0 "0" 0>
  <l l:l1 1 260 200 -26 10 0 0 "" 0 "1e-9" 1 "1e-08" 0 "0" 0>
  <c c:c1 1 340 200 -26 17 0 0 "" 0 "20e-11" 1 "1e-08" 0 "0" 0>
  <l l:l3 1 470 70 -53 -26 0 3 "" 0 "15e-9" 1 "1e-08" 0 "0" 0>
  <LRT * 1 630 270 0 0 0 0>
  <svhb svhb1 1 60 30 0 36 0 0 "3" 1 "f0" 1 "1" 0 "0" 0 "0" 0 "1" 0 "f1" 1 "3" 1 "0" 0 "0" 0 "0" 0 "0" 0 "0" 0>
  <vsource vsource:vs0 1 50 240 18 -26 0 1 "" 0 "0" 0 ".4" 1 "f0" 1 "0" 0 "0" 0 "0" 0>
  <vsource vsource:vs1 1 50 300 18 -26 0 1 "" 0 "0" 0 ".4" 1 "f1" 1 "0" 0 "0" 0 "0" 0>
  <.top top1 1 70 -50 -10 14 0 0 ".options f0=1.0e9 f1=1.001e9 temp=300 gnuplot=1" 1>
  <.out plot1 1 90 410 -35 14 0 0 "Mode=plot" 0 "term" 1 "111" 1 "" 0 "" 0 "vf" 1 "term" 1 "123" 1 "" 0 "" 0 "vf" 1 "sub mag set logscale x; set data style i" 1 "mag.vds" 1>
  <.out plot2 1 220 420 -35 14 0 0 "Mode=plot" 0 "term" 1 "1000" 1 "" 0 "" 0 "vf" 1 "" 0 "" 0 "" 0 "" 0 "" 0 "mag set logscale x; set data style i" 1 "mag.temp" 1>
  <.out plot3 1 350 370 -35 14 0 0 "Mode=plot" 0 "element" 1 "" 0 "mesfetct:m1" 1 "1" 1 "if" 1 "" 0 "" 0 "" 0 "" 0 "" 0 "mag set logscale x; set data style i" 1 "mag.ids" 1>
  <.out plot4 1 620 350 -35 14 0 0 "Mode=plot" 0 "element" 1 "" 0 "mesfetct:m1" 1 "2" 1 "if" 1 "" 0 "" 0 "" 0 "" 0 "" 0 "mag set logscale x; set data style i" 1 "mag.power" 1>
</Components>
<Wires>
  <420 230 440 230 "" 0 0 0 "">
  <420 200 420 230 "" 0 0 0 "">
  <470 260 470 280 "123" 400 250 10 "">
  <260 340 470 340 "" 0 0 0 "">
  <210 200 230 200 "" 0 0 0 "">
  <290 200 310 200 "" 0 0 0 "">
  <370 200 420 200 "" 0 0 0 "">
  <200 60 230 60 "" 0 0 0 "">
  <200 60 200 80 "" 0 0 0 "">
  <350 60 370 60 "" 0 0 0 "">
  <370 60 370 200 "" 0 0 0 "">
  <470 100 470 200 "111" 410 110 41 "">
  <470 100 510 100 "" 0 0 0 "">
  <570 100 590 100 "" 0 0 0 "">
  <650 100 670 100 "" 0 0 0 "">
  <470 -40 530 -40 "" 0 0 0 "">
  <470 -40 470 -20 "" 0 0 0 "">
  <530 -40 530 -10 "" 0 0 0 "">
  <770 250 770 270 "tref" 800 270 9 "">
  <480 180 600 180 "1000" 530 150 19 "">
  <480 180 480 210 "" 0 0 0 "">
  <660 180 770 180 "" 0 0 0 "">
  <770 180 770 190 "" 0 0 0 "">
  <480 270 630 270 "" 0 0 0 "">
  <480 250 480 270 "" 0 0 0 "">
  <630 270 770 270 "" 0 0 0 "">
  <50 200 50 210 "" 0 0 0 "">
  <50 200 150 200 "" 0 0 0 "">
  <50 340 260 340 "" 0 0 0 "">
  <50 330 50 340 "" 0 0 0 "">
</Wires>
<Diagrams>
</Diagrams>
<Paintings>
</Paintings>
