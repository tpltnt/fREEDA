Acknowledgements: 
The gui packaged with fREEDA is based on a gui used in the software simulator Qucs. The Qucs project can be found at http://qucs.sourceforge.net/
The GUI is based on the Qt toolkit (http://qt.nokia.com)

HOW TO ADD AN ELEMENT TO THE GUI
--------------------------------

1. First create the element and properly add it to simulator. Follow the instructions in the README to add a element to the simulator.

2. Create a png representation of your element. Examples of png files can be found in the freeda-<version-number>/simulator/bitmaps directory. The png should be 32x32 pixels.

The next steps will require modifying files in the gui directory.   

3. In the GUI directory you will see a series of files with the prefix fREEDA_. These files represent the categories of elements that will be displayed in the GUI. Decide which category your element belongs to then
open the appropriate file and add the name of your element to arrary.

  For example you have created element Widget and it is a lumped element.
  The Iteminfo of the element in your element class file looks like this
  // Element information
  ItemInfo Resistor::einfo = 
  {
     "wid",
     "Widget",
     "John Smith",
      DEFAULT_ADDRESS"elements/Widget.h.html",
      "2006_07_15"
  };

  Now you would open fREEDA_lumpped.cc which contains the array
  { "abmmixer", "c", "cir", "gyr", "i", "iso", "r", "0"};
  Add the NAME(not neccesarily the class name) but the name of your element.
  The first string in the ItemInfo structure is the name of your elemente,
  which in this example is wid. So now fREEDA_lumped.cc looks like this.
  { "abmmixer", "c", "cir", "gyr", "i", "iso", "r", "wid", "0"};
  Save the change to the file.


4. Open get_info.cpp in the gui directory. Add an if statement exactly like the ones currently in the file however modify the if statement to contain the name of your element.

  So for our widget element we would add these 2 lines
  if (elem_type == Widget::getNetlistName())
  new_elem = new Widget("dummy");


5. Open the interface.cpp file in the gui directory. Follow the instructions in this file to modify it properly.
