/***************************************************************************
                          interface.cpp  -  description
                             -------------------
This file interfaces fREEDA elements and analysis to the gui. If you are adding
an element or analysis to fREEDA and wish to display it in the gui you have to
edit this file. See the drawGenericElement function below.
****************************************************************************/

#include "interface.h"

#include "../simulator/network/NetListItem.h"
#include "../simulator/network/Element.h"
#include "../simulator/network/Circuit.h"

#include "../simulator/analysis/Analysis.h"
#include "../simulator/analysis/Nox_Interface.h"

#include "../simulator/elements/element_headers.h"

#include "../simulator/analysis/AC.h"

#include "../simulator/freeda.h"
#include "../simulator/analysis/anMain.cc"

Interface::Interface(QString& elemtype)
{       
  Element* new_elem = NULL;
  Analysis* new_analysis = NULL; 
  
 //elem_type = elemtype.data();
 // Added by Shivam
   elem_type = elemtype;
  if ((elem_type == "nport4")  || (elem_type == "nport6")  ||
      (elem_type == "nport8")  || (elem_type == "nport10"))
  {     
     nport_str = elem_type;
     elem_type = "nport";
  }
  
 //File needs to be automatically generated  
#include "get_info.cpp"

  if (new_elem != NULL)
  { 
     numberOfTerminals = new_elem->getNumTerms();
     drawGenericElement();
    
     Description = QString(new_elem->getDescription().c_str());      
     Model = new_elem->getName().c_str();
     Name  = Model+":"+Model[0]; //new_elem->getName().c_str();

     numOfPrms = new_elem->getNumberOfParams(); //numparms;
     for (unsigned i = 0; i < numOfPrms; i++)
     {//Append properties to display box
	new_elem->getParamDesc(i, parname, comment, type, dflt_val, required);
	if (dflt_val == "true")
	{//fREEDA parser doesnt like "true" and "false" strings, change to "1" and "0"
	   dflt_val = "1";
	}
	else if (dflt_val == "false")
	{
	  dflt_val = "0";
	}

	if (i == 0)
	{//Add a model name parameter to all elements.
	   if(Model != "nport")
	   {
	      Props.append(new Property("model" , "", false, "Name of the model that will be used for this element"));
	   }
	}  
		   
	//Change required to a bool, what append expects
	if (required == "yes")
	{
           Props.append(new Property(parname.c_str(), "0", true, comment.c_str())); 
	}
	else
	{
           Props.append(new Property(parname.c_str(), dflt_val.c_str(), false, comment.c_str()));
	}
     }
  }
  else if (new_analysis != NULL)
  {
     Description = QString(new_analysis->getName().c_str());     
     Description.append(" analysis");
     Model = new_analysis->getName().c_str();
     Name  = Model;

     drawGenericAnalysis();
     numOfPrms = new_analysis->getNumberOfParams(); //numparms;
     for (unsigned i = 0; i < numOfPrms; i++)
     {//Append properties to display box
	new_analysis->getParamDesc(i, parname, comment, type, dflt_val, required);
        if (dflt_val == "true")
	{//fREEDA parser doesnt like "true" and "false" strings, change to "1" and "0"
	   dflt_val = "1";
	}
	else if (dflt_val == "false")
	{
	  dflt_val = "0";
	}

	//Change required to a bool, what append expects
	if (required == "yes")
	{
	   Props.append(new Property(parname.c_str(), "0", true, comment.c_str())); 
	}
	else
	{
	   Props.append(new Property(parname.c_str(), dflt_val.c_str(), false, comment.c_str()));
	}	
     } 
  }
  else
     return;
  		
}

Interface::~Interface()
{
}

void Interface::drawGenericAnalysis()
{
  QString s = Description;
  int a = s.find(" ");
  int b = s.findRev(" ");
  if (a != -1 && b != -1) {
    if (a > (int) s.length() - b)  b = a;
  }
  if (a < 8 || s.length() - b < 8) b = -1;
  if (b != -1) s[b] = '\n';

  Lines.append(new Line(  0,  0,100,  0,QPen(Qt::darkBlue,2)));
  Lines.append(new Line(  0, 40,100, 40,QPen(Qt::darkBlue,2)));
  Lines.append(new Line(  0,  0,  0, 40,QPen(Qt::darkBlue,2)));
  Lines.append(new Line(100,  0,100, 40,QPen(Qt::darkBlue,2)));

  Texts.append(new Text(0, 0, s.left(b)));
  if (b != -1) Texts.append(new Text(0, 0, s.mid(b+1)));

  x1 = -10; y1 = -9;
  x2 = x1+128; y2 = y1+41;

  tx = 0;
  ty = y2+1;
}

/***************************************************************
The following fucntion will draw the element selected on to the
background of the gui. If you are adding an element you need
to add a else if statement to this function. For example I am adding
element widget, then below in the function I would add.

  else if (elem_type == "widget") 
  {  
    //These commands draw your element. You can draw lines and circles.
    //See the README for more detail
    Lines.append(new Line(-10,-15,-10, 15,QPen(Qt::darkBlue,2)));
    Lines.append(new Line(-30,  0,-10,  0,QPen(Qt::darkBlue,2)));
    Lines.append(new Line(-10, -5,  0,-15,QPen(Qt::darkBlue,2)));
    Arcs.append(new Arc(-14,-14, 29, 29,  0,16*360,QPen(Qt::darkBlue,2)));
    

    //These commands give the locations of the ports into your element.
    Ports.append(new Port(-30,  0));
    Ports.append(new Port( 30,  0));

    //Give the corners of your element and specify the location of any
    //text you put to describe you element.    
    x1 = -30; y1 = -16;
    x2 =  30; y2 =  30;

    tx = x1+4;
    ty = y1-5;
  }

For documentation on on the Lines, Arcs and Ports classes please see
element_qucs.h and Qt documentation and the QPainter class. 


****************************************************************/
void Interface::drawGenericElement()
{

  if ((elem_type == "abmvcann")  || (elem_type == "abmvtanh"))
  {
    Arcs.append(new Arc(0,-18, 23, 23,  0, 16*360,QPen(Qt::darkBlue,2)));
    Lines.append(new Line( 11,-14, 11,  0,QPen(Qt::darkBlue,3)));
    Lines.append(new Line( 11, -1, 15, -6,QPen(Qt::darkBlue,3)));
    Lines.append(new Line( 11, -1,  7, -6,QPen(Qt::darkBlue,3)));

    Lines.append(new Line(-30,-30,-12,-30,QPen(Qt::darkBlue,2)));
    Lines.append(new Line( 11,-30, 30,-30,QPen(Qt::darkBlue,2)));
    Lines.append(new Line( 11, 16,-12, 16,QPen(Qt::darkBlue,2)));
    Lines.append(new Line(  0, 16,  0, 30,QPen(Qt::darkBlue,2)));

    Lines.append(new Line(-12,-30,-12,-23,QPen(Qt::darkBlue,2)));
    Lines.append(new Line(-12, 16,-12, 9,QPen(Qt::darkBlue,2)));
    Lines.append(new Line( 11,-30, 11,-18,QPen(Qt::darkBlue,2)));
    Lines.append(new Line( 11, 16, 11, 4,QPen(Qt::darkBlue,2)));

    Arcs.append(new Arc(-16,-23, 8, 8,  0, 16*360,QPen(Qt::darkBlue,2)));
    Arcs.append(new Arc(-16,  1, 8, 8,  0, 16*360,QPen(Qt::darkBlue,2)));
    
    Lines.append(new Line(-25,-27, 25,-27,QPen(Qt::darkGray,1)));
    Lines.append(new Line( 25,-27, 25, 24,QPen(Qt::darkGray,1)));
    Lines.append(new Line( 25, 24,-25, 24,QPen(Qt::darkGray,1)));
    Lines.append(new Line(-25, 24,-25,-27,QPen(Qt::darkGray,1)));

    Lines.append(new Line( 19,-24, 19,-18,QPen(Qt::red,1)));
    Lines.append(new Line( 16,-21, 22,-21,QPen(Qt::red,1)));
    
    Lines.append(new Line( 16, 10, 22, 10,QPen(Qt::black,1)));

    Ports.append(new Port(-30,-30));
    Ports.append(new Port( 30,-30));
    Ports.append(new Port( 0, 30));


    x1 = -30; y1 = -30;
    x2 =  30; y2 =  30;

    tx = x1+4;
    ty = y2+4;
  
  }
  else if (elem_type == "abmmixer")
  {//Variable inputs
    Arcs.append(new Arc(-14,-14, 29, 29,  0,16*360,QPen(Qt::darkBlue,2)));
    Lines.append(new Line(-30,  0,-14,  0,QPen(Qt::darkBlue,2)));
    Lines.append(new Line( 30,  0, 14,  0,QPen(Qt::darkBlue,2)));
    Lines.append(new Line(  0, 14,  0, 30,QPen(Qt::darkBlue,2)));

    Lines.append(new Line(-8, -10, 8, 10,QPen(Qt::darkBlue,2)));
    Lines.append(new Line(-8, 10, 8, -10,QPen(Qt::darkBlue,2)));
       
    Lines.append(new Line(-22, -4, -26,  4,QPen(Qt::darkBlue,2)));   // marks port 1

    Ports.append(new Port(-30,  0));
    Ports.append(new Port( 30,  0));
    Ports.append(new Port(  0, 30));

    x1 = -30; y1 = -16;
    x2 =  30; y2 =  30;

    tx = x1+4;
    ty = y2+5;
  }
  else if ((elem_type == "abmbutterbpf10") || (elem_type == "attenuatorgen") || (elem_type == "zbutter") || (elem_type == "zchebyi"))   
  {
    Lines.append(new Line(-14,-14, 14,-14,QPen(Qt::darkBlue,2)));
    Lines.append(new Line(-14, 14, 14, 14,QPen(Qt::darkBlue,2)));
    Lines.append(new Line(-14,-14,-14, 14,QPen(Qt::darkBlue,2)));
    Lines.append(new Line( 14,-14, 14, 14,QPen(Qt::darkBlue,2)));

    Lines.append(new Line(-30, 10,-14, 10,QPen(Qt::darkBlue,2)));
    Lines.append(new Line(-30,-10,-14,-10,QPen(Qt::darkBlue,2)));    
    Lines.append(new Line( 14,-10, 30,-10,QPen(Qt::darkBlue,2)));
    Lines.append(new Line( 14, 10, 30, 10,QPen(Qt::darkBlue,2)));

    Ports.append(new Port(-30,-10));
    Ports.append(new Port(-30, 10));
    Ports.append(new Port( 30,-10));
    Ports.append(new Port( 30, 10));
    
    x1 = -30; y1 = -17;
    x2 =  30; y2 =  17;

    tx = x1+4;
    ty = y2+4;
  
  } 
  else if ((elem_type == "bjtnpn")  || 
           (elem_type == "bjtnpnt") || 
	   (elem_type == "bjtpnp")  ||
	   (elem_type == "hbtnpnxt"))
  {
  
    Lines.append(new Line(-10,-15,-10, 15,QPen(Qt::darkBlue,3)));
    Lines.append(new Line(-30,  0,-10,  0,QPen(Qt::darkBlue,2)));
    Lines.append(new Line(-10, -5,  0,-15,QPen(Qt::darkBlue,2)));
    Lines.append(new Line(  0,-15,  0,-30,QPen(Qt::darkBlue,2)));
    Lines.append(new Line(-10,  5,  0, 15,QPen(Qt::darkBlue,2)));
    Lines.append(new Line(  0, 15,  0, 30,QPen(Qt::darkBlue,2)));

    Lines.append(new Line(  9,  0, 30,  0,QPen(Qt::darkBlue,2)));
    Lines.append(new Line(  9, -7,  9,  7,QPen(Qt::darkBlue,3)));

    Lines.append(new Line( -6, 15,  0, 15,QPen(Qt::darkBlue,2)));
    Lines.append(new Line(  0,  9,  0, 15,QPen(Qt::darkBlue,2)));

    Ports.append(new Port(  0,-30));
    Ports.append(new Port(-30,  0));
    Ports.append(new Port(  0, 30));
    Ports.append(new Port( 30,  0));
  
    if ((elem_type == "bjtnpnt") || (elem_type == "hbtnpnxt"))
    {
       Ports.append(new Port( 30,-20));
       Ports.append(new Port( 30, 20));
       Texts.append(new Text( 10, 10,"ref"));
    }
  

    x1 = -30; y1 = -30;
    x2 =  30; y2 =  30;

    tx = x2+15;
    ty = y1+4;
    Line *pl2 = Lines.last();
    Line *pl1 = Lines.prev();

    if(elem_type != "bjtpnp") {
      pl1->x1 = -6;  pl1->y1 = 15;  pl1->x2 = 0; pl1->y2 = 15;
      pl2->x1 = 0;  pl2->y1 = 9;  pl2->x2 = 0; pl2->y2 = 15;
    }
    else {
      pl1->x1 = -5;  pl1->y1 = 10;  pl1->x2 = -5; pl1->y2 = 16;
      pl2->x1 = -5;  pl2->y1 = 10;  pl2->x2 = 1; pl2->y2 = 10;
    }
  }
  else if ((elem_type == "butterworthbpf") || (elem_type == "chebyshevbpf") || (elem_type == "chebyshevbsf") || (elem_type == "chebyshevhpf") || (elem_type == "chebyshevlpf") || (elem_type == "lspiral"))
  {
        Lines.append(new Line(-14,-14, 14,-14,QPen(Qt::darkBlue,2)));
        Lines.append(new Line(-14, 14, 14, 14,QPen(Qt::darkBlue,2)));
        Lines.append(new Line(-14,-14,-14, 14,QPen(Qt::darkBlue,2)));
        Lines.append(new Line( 14,-14, 14, 14,QPen(Qt::darkBlue,2)));

        Lines.append(new Line(-30, 0,-14, 0,QPen(Qt::darkBlue,2)));
        Lines.append(new Line( 14, 0, 30, 0,QPen(Qt::darkBlue,2)));
        Lines.append(new Line(  0, 14,  0, 30,QPen(Qt::darkBlue,2)));

        Ports.append(new Port(-30,0));
        Ports.append(new Port(30,0));
        Ports.append(new Port(  0, 30));

            x1 = -30; y1 = -17;
            x2 =  30; y2 =  17;

            tx = x1+4;
            ty = y2+4;
  }

  else if ((elem_type == "c")     ||
           (elem_type == "capacitor")   ||
           (elem_type == "capacitorjn") ||
	   (elem_type == "capacitorlc") ||
	   (elem_type == "capacitormos") || (elem_type == "capacitorferroelectric"))  
  {
    Lines.append(new Line( -4,-11, -4, 11,QPen(Qt::darkBlue,4)));
    Lines.append(new Line(  4,-11,  4, 11,QPen(Qt::darkBlue,4)));
    Lines.append(new Line(-30,  0, -4,  0,QPen(Qt::darkBlue,2)));
    Lines.append(new Line(  4,  0, 30,  0,QPen(Qt::darkBlue,2)));

    Ports.append(new Port(-30,  0));
    Ports.append(new Port( 30,  0));

    x1 = -30; y1 = -13;
    x2 =  30; y2 =  13;

    tx = x1+4;
    ty = y2+4; 
  }
  else if (elem_type == "circulator")
  {
    Arcs.append(new Arc(-14,-14, 29, 29,  0,16*360,QPen(Qt::darkBlue,2)));
    Lines.append(new Line(-30,  0,-14,  0,QPen(Qt::darkBlue,2)));
    Lines.append(new Line( 30,  0, 14,  0,QPen(Qt::darkBlue,2)));
    Lines.append(new Line(  0, 14,  0, 30,QPen(Qt::darkBlue,2)));
    Lines.append(new Line(  0, 14,  0,-30,QPen(Qt::darkBlue,2)));

    Arcs.append(new Arc( -8, -6, 17, 17,16*20,16*150,QPen(Qt::darkBlue,2)));
    Lines.append(new Line(  8,  0,  9, -7,QPen(Qt::darkBlue,2)));
    Lines.append(new Line(  8,  0,  2, -1,QPen(Qt::darkBlue,2)));

    Lines.append(new Line(-22, -4,-26,  4,QPen(Qt::darkBlue,2)));   // marks port 1

    Ports.append(new Port(-30,  0));
    Ports.append(new Port( 30,  0));
    Ports.append(new Port(  0, 30));

    Ports.append(new Port(  0,-30));
    Ports.append(new Port(  0, 30));
    Ports.append(new Port(  0, 30));
 
    x1 = -30; y1 = -16;
    x2 =  30; y2 =  30;

    tx = x1+4;
    ty = y2+5;
  }
  else if (elem_type == "cpw")
  {
  
    Lines.append(new Line(-30,-10,-22,-10,QPen(Qt::darkBlue,2)));
    Lines.append(new Line( 22,-10, 30,-10,QPen(Qt::darkBlue,2)));
    Lines.append(new Line(-30, 10,-22, 10,QPen(Qt::darkBlue,2)));
    Lines.append(new Line( 22, 10, 30, 10,QPen(Qt::darkBlue,2)));
   
    Lines.append(new Line(-22,-18, 22,-18,QPen(Qt::darkBlue,2)));
    Lines.append(new Line(-22,-18,-22, -2,QPen(Qt::darkBlue,2)));
    Lines.append(new Line(-22, -2, 22, -2,QPen(Qt::darkBlue,2)));
    Lines.append(new Line( 22,-18, 22, -2,QPen(Qt::darkBlue,2)));
    
    Lines.append(new Line(-22, 18, 22, 18,QPen(Qt::darkBlue,2)));
    Lines.append(new Line(-22, 18,-22,  2,QPen(Qt::darkBlue,2)));
    Lines.append(new Line(-22,  2, 22,  2,QPen(Qt::darkBlue,2)));
    Lines.append(new Line( 22, 18, 22,  2,QPen(Qt::darkBlue,2)));
    
    Lines.append(new Line(-25,-21, 25,-21,QPen(Qt::darkBlue,2)));
    Lines.append(new Line( 16,-29, 24,-21,QPen(Qt::darkBlue,2)));
    Lines.append(new Line(  8,-29, 16,-21,QPen(Qt::darkBlue,2)));
    Lines.append(new Line(  0,-29,  8,-21,QPen(Qt::darkBlue,2)));
    Lines.append(new Line( -8,-29,  0,-21,QPen(Qt::darkBlue,2)));
    Lines.append(new Line(-16,-29, -8,-21,QPen(Qt::darkBlue,2)));
    Lines.append(new Line(-24,-29,-16,-21,QPen(Qt::darkBlue,2)));

    Lines.append(new Line(-25, 21, 25, 21,QPen(Qt::darkBlue,2)));
    Lines.append(new Line(-24, 21,-16, 29,QPen(Qt::darkBlue,2)));
    Lines.append(new Line(-16, 21, -8, 29,QPen(Qt::darkBlue,2)));
    Lines.append(new Line( -8, 21,  0, 29,QPen(Qt::darkBlue,2)));
    Lines.append(new Line(  0, 21,  8, 29,QPen(Qt::darkBlue,2)));
    Lines.append(new Line(  8, 21, 16, 29,QPen(Qt::darkBlue,2)));
    Lines.append(new Line( 16, 21, 24, 29,QPen(Qt::darkBlue,2)));

    Ports.append(new Port(-30,-10));
    Ports.append(new Port(-30, 10));
    Ports.append(new Port( 30,-10));
    Ports.append(new Port( 30, 10));

    x1 = -30; y1 =-24;
    x2 =  30; y2 = 24;

    tx = x1+4;
    ty = y2+4;
  }
  else if (elem_type == "dhld")
  {
//    Lines.append(new Line( 11, -7, 11,  7,QPen(Qt::darkBlue,3)));
//    Lines.append(new Line( 11,  6, 15,  1,QPen(Qt::darkBlue,3)));
//    Lines.append(new Line( 11,  6,  7,  1,QPen(Qt::darkBlue,3)));

    Lines.append(new Line(-30,-20, -6,-20,QPen(Qt::darkBlue,2)));
    Lines.append(new Line(-30, 20, -7, 20,QPen(Qt::darkBlue,2)));
    Lines.append(new Line( 30,-20, 20,-20,QPen(Qt::darkBlue,2)));
    Lines.append(new Line( 30, 20, 20, 20,QPen(Qt::darkBlue,2)));

//    Lines.append(new Line(-12,-30,-12,-23,QPen(Qt::darkBlue,2)));
    Lines.append(new Line( -6,-20, -6,-10,QPen(Qt::darkBlue,2)));
    Lines.append(new Line( -7, 20, -7,  7,QPen(Qt::darkBlue,2)));

    Lines.append(new Line(-12,-10, -1,-10,QPen(Qt::darkBlue,2)));
    Lines.append(new Line(-12,-10, -7,  7,QPen(Qt::darkBlue,2)));
    Lines.append(new Line( -1,-10, -7,  7,QPen(Qt::darkBlue,2)));
    Lines.append(new Line(-10,  7, -4,  7,QPen(Qt::darkBlue,2)));

    Lines.append(new Line(-20,-27, 20,-27,QPen(Qt::darkBlue,2)));
    Lines.append(new Line( 20,-27, 20, 27,QPen(Qt::darkBlue,2)));
    Lines.append(new Line( 20, 27,-20, 27,QPen(Qt::darkBlue,2)));
    Lines.append(new Line(-20, 27,-20,-27,QPen(Qt::darkBlue,2)));

    Ports.append(new Port(-30,-20));
    Ports.append(new Port(-30, 20));
    Ports.append(new Port( 30,-20));
    Ports.append(new Port( 30, 20));


    x1 = -30; y1 = -30;
    x2 =  30; y2 =  30;

    tx = x1+4;
    ty = y2+4;
  }  
  
  else if ((elem_type == "d")            || 
           (elem_type == "diode")        || 
	   (elem_type == "diodec")       || 
	   (elem_type == "diodecompact") ||
     (elem_type == "diodeqk")      ||
	   (elem_type == "diodesp")      ||	   
	   (elem_type == "diodetun"))
  {
    Lines.append(new Line(-30,  0, -6,  0,QPen(Qt::darkBlue,2)));
    Lines.append(new Line( 30,  0,  6,  0,QPen(Qt::darkBlue,2)));
    Lines.append(new Line( -6, -9, -6,  9,QPen(Qt::darkBlue,2)));
    Lines.append(new Line(  6, -9,  6,  9,QPen(Qt::darkBlue,2)));
    Lines.append(new Line( -6, -9,  6,  0,QPen(Qt::darkBlue,2)));
    Lines.append(new Line( -6,  9,  6,  0,QPen(Qt::darkBlue,2)));

    Ports.append(new Port(-30, 0));
    Ports.append(new Port( 30, 0));
    
    x1 = -30; y1 = -9;
    x2 =  30; y2 =  9;

    tx = x1+4;
    ty = y2+4;
  }
  
  else if((elem_type == "dflipflop") || (elem_type == "dflipflopx"))
  {
        Lines.append(new Line(-14,-16,-14, 16,QPen(Qt::darkBlue,2)));
        Lines.append(new Line(-14,-16, 14 ,-16,QPen(Qt::darkBlue,2)));
        Lines.append(new Line(-14, 16, 14, 16,QPen(Qt::darkBlue,2)));
        Lines.append(new Line(14, -16, 14, 16,QPen(Qt::darkBlue,2)));

        Arcs.append(new Arc( 14, 6,  7,  7,  0,16*360,QPen(Qt::darkBlue,2)));
        Arcs.append(new Arc( -21, 6,  7,  7,  0,16*360,QPen(Qt::darkBlue,2)));

        Lines.append(new Line(-30, 10,-21, 10,QPen(Qt::darkBlue,2)));
        Lines.append(new Line(-30,-10,-14,-10,QPen(Qt::darkBlue,2)));
        Lines.append(new Line(-14,12,-6,9,QPen(Qt::darkBlue,2)));
        Lines.append(new Line(-14,6,-6,9,QPen(Qt::darkBlue,2)));
        Lines.append(new Line(21,10,30,10,QPen(Qt::darkBlue,2)));
        Lines.append(new Line(14,-10,30,-10,QPen(Qt::darkBlue,2)));

        Lines.append(new Line(  0, 16,  0, 30,QPen(Qt::darkBlue,2)));
        Lines.append(new Line(  0,-16,  0,-30,QPen(Qt::darkBlue,2)));

        Ports.append(new Port(  0,-30));
        Ports.append(new Port(-30,-10));
        Ports.append(new Port(-30, 10));
        Ports.append(new Port( 30, 10));
        Ports.append(new Port( 30, -10));
        Ports.append(new Port(  0, 30));

        x1 = -30; y1 = -30;
        x2 =  30; y2 =  30;

        tx = x2-8;
        ty = y2;
  }


  else if (elem_type == "gyrator")
  {
     Arcs.append(new Arc(  2, -9, 19, 19, 16*90, 16*180,QPen(Qt::darkBlue,2)));
     Arcs.append(new Arc(-21, -9, 19, 19,16*270, 16*180,QPen(Qt::darkBlue,2)));

     Lines.append(new Line(-30,-30,-12,-30,QPen(Qt::darkBlue,2)));
     Lines.append(new Line(-30, 30,-12, 30,QPen(Qt::darkBlue,2)));
     Lines.append(new Line( 12,-30, 30,-30,QPen(Qt::darkBlue,2)));
     Lines.append(new Line( 12, 30, 30, 30,QPen(Qt::darkBlue,2)));

     Lines.append(new Line( 12,-30, 12, 30,QPen(Qt::darkBlue,2)));
     Lines.append(new Line(-12,-30,-12, 30,QPen(Qt::darkBlue,2)));

     Lines.append(new Line(-22,-22, 22,-22,QPen(Qt::darkGray,1)));
     Lines.append(new Line( 22,-22, 22, 22,QPen(Qt::darkGray,1)));
     Lines.append(new Line( 22, 22,-22, 22,QPen(Qt::darkGray,1)));
     Lines.append(new Line(-22, 22,-22,-22,QPen(Qt::darkGray,1)));

     Ports.append(new Port(-30,-30));
     Ports.append(new Port(-30, 30));     
     Ports.append(new Port( 30,-30));
     Ports.append(new Port( 30, 30));

     x1 = -30; y1 = -30;
     x2 =  30; y2 =  30;

     tx = x1+4;
     ty = y2+4;  
  }
  else if (elem_type == "isolator")
  {
    Lines.append(new Line( -8,  0,  8,  0,QPen(Qt::darkBlue,3)));
    Lines.append(new Line(  8,  0,  0, -5,QPen(Qt::darkBlue,3)));
    Lines.append(new Line(  8,  0,  0,  5,QPen(Qt::darkBlue,3)));

    Lines.append(new Line(-14,-14, 14,-14,QPen(Qt::darkBlue,2)));
    Lines.append(new Line(-14, 14, 14, 14,QPen(Qt::darkBlue,2)));
    Lines.append(new Line(-14,-14,-14, 14,QPen(Qt::darkBlue,2)));
    Lines.append(new Line( 14,-14, 14, 14,QPen(Qt::darkBlue,2)));

    Lines.append(new Line(-30,-10,-14,-10,QPen(Qt::darkBlue,2)));
    Lines.append(new Line( 14,-10, 30,-10,QPen(Qt::darkBlue,2)));
    Lines.append(new Line(-30, 10,-14, 10,QPen(Qt::darkBlue,2)));
    Lines.append(new Line( 14, 10, 30, 10,QPen(Qt::darkBlue,2)));

    Ports.append(new Port(-30, -10));
    Ports.append(new Port( 30, -10));
    Ports.append(new Port(-30,  10));
    Ports.append(new Port( 30,  10));

    x1 = -30; y1 = -17;
    x2 =  30; y2 =  17;

    tx = x1+4;
    ty = y2+4;
  }
  else if ((elem_type == "inductor") || 
           (elem_type == "i")        || 
	   (elem_type == "l"))
  {
    Arcs.append(new Arc(-18, -6, 13, 13,  0, 16*180,QPen(Qt::darkBlue,2)));
    Arcs.append(new Arc( -6, -6, 13, 13,  0, 16*180,QPen(Qt::darkBlue,2)));
    Arcs.append(new Arc(  6, -6, 13, 13,  0, 16*180,QPen(Qt::darkBlue,2)));
    Lines.append(new Line(-30,  0,-18,  0,QPen(Qt::darkBlue,2)));
    Lines.append(new Line( 18,  0, 30,  0,QPen(Qt::darkBlue,2)));

    Ports.append(new Port(-30,  0));
    Ports.append(new Port( 30,  0));

    x1 = -30; y1 = -10;
    x2 =  30; y2 =   6;

    tx = x1+4;
    ty = y2+4;
  }

  else if ((elem_type == "interconnectrt") || (elem_type == "interconnectrtsh") || (elem_type == "viart")) 
  {
        Lines.append(new Line(-30,  0,-18,  0,QPen(Qt::darkBlue,2)));
        Lines.append(new Line(-18,  0,-15, -7,QPen(Qt::darkBlue,2)));
        Lines.append(new Line(-15, -7, -9,  7,QPen(Qt::darkBlue,2)));
        Lines.append(new Line( -9,  7, -3, -7,QPen(Qt::darkBlue,2)));
        Lines.append(new Line( -3, -7,  3,  7,QPen(Qt::darkBlue,2)));
        Lines.append(new Line(  3,  7,  9, -7,QPen(Qt::darkBlue,2)));
        Lines.append(new Line(  9, -7, 15,  7,QPen(Qt::darkBlue,2)));
        Lines.append(new Line( 15,  7, 18,  0,QPen(Qt::darkBlue,2)));
        Lines.append(new Line( 18,  0, 30,  0,QPen(Qt::darkBlue,2)));

        Lines.append(new Line(-18,-7, 18,-7,QPen(Qt::darkGray,1)));
        Lines.append(new Line( 18,-7, 18, 7,QPen(Qt::darkGray,1)));
        Lines.append(new Line( 18, 7,-18, 7,QPen(Qt::darkGray,1)));
        Lines.append(new Line(-18, 7,-18,-7,QPen(Qt::darkGray,1)));

        Ports.append(new Port(-30,  0));
        Ports.append(new Port( 30,  0));
        Ports.append(new Port( 15, 15));
        Texts.append(new Text( 17, 14,"Tref"));
        Ports.append(new Port(-15, 15));

        x1 = -30; y1 = -11;
        x2 =  30; y2 =  22;

        tx = x1+4;
        ty = y2+4;

  }

  else if ((elem_type == "isffm") ||
           (elem_type == "iam")      ||
	   (elem_type == "iexp")     ||
	   (elem_type == "isource"))
  {
    Arcs.append(new Arc(-12,-12, 25, 25,  0, 16*360,QPen(Qt::darkBlue,2)));
    Lines.append(new Line(-30,  0,-12,  0,QPen(Qt::darkBlue,2)));
    Lines.append(new Line( 30,  0, 12,  0,QPen(Qt::darkBlue,2)));
    Lines.append(new Line( -7,  0,  7,  0,QPen(Qt::darkBlue,3)));
    Lines.append(new Line(  6,  0,  0, -4,QPen(Qt::darkBlue,3)));
    Lines.append(new Line(  6,  0,  0,  4,QPen(Qt::darkBlue,3)));
    Arcs.append(new Arc( 12,  5,  7,  7,16*270, 16*180,QPen(Qt::darkBlue,2)));
    Arcs.append(new Arc( 12, 11,  7,  7, 16*90, 16*180,QPen(Qt::darkBlue,2)));

    Ports.append(new Port( 30,  0));
    Ports.append(new Port(-30,  0));

    x1 = -30; y1 = -14;
    x2 =  30; y2 =  16;

    tx = x1+4;
    ty = y2+4;
  }
  else if ((elem_type == "jfetn")    || (elem_type == "jfetp")    ||
           (elem_type == "mesfetc")  || (elem_type == "mesfetcq") ||
	   (elem_type == "mesfetct") || (elem_type == "mesfetm")  ||
	   (elem_type == "mesfetmq") || (elem_type == "mesfettom"))
  {
    Lines.append(new Line(-10,-15,-10, 15,QPen(Qt::darkBlue,3)));
    Lines.append(new Line(-30,  0,-10,  0,QPen(Qt::darkBlue,2)));
    Lines.append(new Line(-10,-10,  0,-10,QPen(Qt::darkBlue,2)));
    Lines.append(new Line(  0,-10,  0,-30,QPen(Qt::darkBlue,2)));
    Lines.append(new Line(-10, 10,  0, 10,QPen(Qt::darkBlue,2)));
    Lines.append(new Line(  0, 10,  0, 30,QPen(Qt::darkBlue,2)));

    Lines.append(new Line( -4, 24,  4, 20,QPen(Qt::darkBlue,2)));

    // these two lines must be the last ones !
    Lines.append(new Line(-16, -5,-11,  0,QPen(Qt::darkBlue,2)));
    Lines.append(new Line(-16,  5,-11,  0,QPen(Qt::darkBlue,2)));

    if ((elem_type == "jfetn")  || (elem_type == "jfetp") || 
        (elem_type == "mosntft"))
    {
       Ports.append(new Port(  0,-30));
       Ports.append(new Port(-30,  0));
       Ports.append(new Port(  0, 30));
    }
    if ((elem_type == "mesfetc")  || (elem_type == "mesfetcq") ||
	(elem_type == "mesfetct") || (elem_type == "mesfetm")  ||
	(elem_type == "mesfetmq") || (elem_type == "mesfettom"))
    {
       Ports.append(new Port(-30,  0));
       Ports.append(new Port(  0,-30));
       Ports.append(new Port(  0, 30));    
    }
    if (elem_type == "mesfetct")
    {
       Ports.append(new Port( 10,-20));
       Ports.append(new Port( 10, 20));
       Texts.append(new Text( 15, 15,"ref"));    
    }


    x1 = -30; y1 = -30;
    x2 =  10; y2 =  30;

    tx = x2+20;
    ty = y1+4;
    
    Line *pl2 = Lines.last();
    Line *pl1 = Lines.prev();

    if(elem_type != "jfetp") {
      pl1->x1 = -16;  pl1->y1 = -5;  pl1->x2 = -11; pl1->y2 = 0;
      pl2->x1 = -16;  pl2->y1 =  5;  pl2->x2 = -11; pl2->y2 = 0;
    }
    else {
      pl1->x1 = -18;  pl1->y1 = 0;  pl1->x2 = -13; pl1->y2 = -5;
      pl2->x1 = -18;  pl2->y1 = 0;  pl2->x2 = -13; pl2->y2 =  5;
    }

  } 
  else if ((elem_type == "k") || (elem_type == "transformer"))
  {
    Arcs.append(new Arc(-16,-18,13,13, 16*270,16*180, QPen(Qt::darkBlue,2)));
    Arcs.append(new Arc(-16, -6,13,13, 16*270,16*180, QPen(Qt::darkBlue,2)));
    Arcs.append(new Arc(-16,  6,13,13, 16*270,16*180, QPen(Qt::darkBlue,2)));
    Arcs.append(new Arc(  4,-18,13,13,  16*90,16*180, QPen(Qt::darkBlue,2)));
    Arcs.append(new Arc(  4, -6,13,13,  16*90,16*180, QPen(Qt::darkBlue,2)));
    Arcs.append(new Arc(  4,  6,13,13,  16*90,16*180, QPen(Qt::darkBlue,2)));
    Lines.append(new Line(-10,-18,-10,-30,QPen(Qt::darkBlue,2)));
    Lines.append(new Line(-10,-30,-30,-30,QPen(Qt::darkBlue,2)));
    Lines.append(new Line( 10,-18, 10,-30,QPen(Qt::darkBlue,2)));
    Lines.append(new Line( 10,-30, 30,-30,QPen(Qt::darkBlue,2)));
    Lines.append(new Line(-10, 18,-10, 30,QPen(Qt::darkBlue,2)));
    Lines.append(new Line(-10, 30,-30, 30,QPen(Qt::darkBlue,2)));
    Lines.append(new Line( 10, 18, 10, 30,QPen(Qt::darkBlue,2)));
    Lines.append(new Line( 10, 30, 30, 30,QPen(Qt::darkBlue,2)));
    Lines.append(new Line( -1,-20, -1, 20,QPen(Qt::darkBlue,1)));
    Lines.append(new Line(  1,-20,  1, 20,QPen(Qt::darkBlue,1)));

    Texts.append(new Text(-21, -18,"T"));
    Arcs.append(new Arc(-21,-24,  6,  6,  0, 16*360,QPen(Qt::darkBlue,2)));
    Arcs.append(new Arc( 15,-24,  6,  6,  0, 16*360,QPen(Qt::darkBlue,2)));


    Ports.append(new Port(-30,-30));
    Ports.append(new Port(-30, 30));
    Ports.append(new Port( 30,-30));
    Ports.append(new Port( 30, 30));

    x1 = -33; y1 = -33;
    x2 =  33; y2 =  33;

    tx = x1+4;
    ty = y2+4;
  }
  else if (elem_type == "transformerct")
  {
    Arcs.append(new Arc(-16,-18,13,13, 16*270,16*180, QPen(Qt::darkBlue,2)));
    Arcs.append(new Arc(-16, -6,13,13, 16*270,16*180, QPen(Qt::darkBlue,2)));
    Arcs.append(new Arc(-16,  6,13,13, 16*270,16*180, QPen(Qt::darkBlue,2)));
    Arcs.append(new Arc(  4,-18,13,13,  16*90,16*180, QPen(Qt::darkBlue,2)));
    Arcs.append(new Arc(  4, -6,13,13,  16*90,16*180, QPen(Qt::darkBlue,2)));
    Arcs.append(new Arc(  4,  6,13,13,  16*90,16*180, QPen(Qt::darkBlue,2)));
    Lines.append(new Line(-10,-18,-10,-30,QPen(Qt::darkBlue,2)));
    Lines.append(new Line(-10,-30,-30,-30,QPen(Qt::darkBlue,2)));
    Lines.append(new Line( 10,-18, 10,-30,QPen(Qt::darkBlue,2)));
    Lines.append(new Line( 10,-30, 30,-30,QPen(Qt::darkBlue,2)));
    Lines.append(new Line(-10, 18,-10, 30,QPen(Qt::darkBlue,2)));
    Lines.append(new Line(-10, 30,-30, 30,QPen(Qt::darkBlue,2)));
    Lines.append(new Line( 10, 18, 10, 30,QPen(Qt::darkBlue,2)));
    Lines.append(new Line( 10, 30, 30, 30,QPen(Qt::darkBlue,2)));
    Lines.append(new Line( -1,-20, -1, 20,QPen(Qt::darkBlue,1)));
    Lines.append(new Line(  1,-20,  1, 20,QPen(Qt::darkBlue,1)));
    Lines.append(new Line( 7,0, 30,0,QPen(Qt::darkBlue,2)));

   Texts.append(new Text(-21, -18,"T"));
   Arcs.append(new Arc(-21,-24,  6,  6,  0, 16*360,QPen(Qt::darkBlue,2)));
   Arcs.append(new Arc( 15,-24,  6,  6,  0, 16*360,QPen(Qt::darkBlue,2)));
   Arcs.append(new Arc( 15,15,  6,  6,  0, 16*360,QPen(Qt::darkBlue,2)));

   Ports.append(new Port(-30,-30));
   Ports.append(new Port(-30, 30));
   Ports.append(new Port( 30,-30));
   Ports.append(new Port( 30, 30));
   Ports.append(new Port( 30, 0));

       x1 = -33; y1 = -33;
       x2 =  33; y2 =  33;

           tx = x1+4;
           ty = y2+4;
}


  else if ((elem_type == "molecule1") || (elem_type == "molecule2"))
  {
    Arcs.append(new Arc(-14,-14, 25, 25,  0,16*360,QPen(Qt::darkBlue,2)));
    Lines.append(new Line(-36,  0,-20,  0,QPen(Qt::darkBlue,2)));
    Lines.append(new Line( 33,  0, 17,  0,QPen(Qt::darkBlue,2)));

    Arcs.append(new Arc(-20, -3, 6, 6,  0,16*360,QPen(Qt::darkBlue,2)));
    Arcs.append(new Arc( 11, -3, 6, 6,  0,16*360,QPen(Qt::darkBlue,2)));
    Lines.append(new Line(-24, -4,-28,  4,QPen(Qt::darkBlue,2)));   // marks port 1

    Ports.append(new Port(-30,  0));
    Ports.append(new Port( 30,  0));

    x1 = -30; y1 = -16;
    x2 =  30; y2 =  16;

    tx = x1+4;
    ty = y2+4;  
  }
  else if ((elem_type == "cmosinv")   ||
  	   (elem_type == "cmosinvt")  ||
	   (elem_type == "opampt"))
  {
    Lines.append(new Line(-14,-20,-14, 20,QPen(Qt::darkBlue,2)));
    Lines.append(new Line(-14,-20, 14,  0,QPen(Qt::darkBlue,2)));
    Lines.append(new Line(-14, 20, 14,  0,QPen(Qt::darkBlue,2)));
    Arcs.append(new Arc( 14, -3, 6, 6,  0,16*360,QPen(Qt::darkBlue,2)));
    
    Lines.append(new Line(-30,  0,-14,  0,QPen(Qt::darkBlue,2)));
    Lines.append(new Line( 20,  0, 30,  0,QPen(Qt::darkBlue,2)));
    
    Lines.append(new Line(  0, 10,  0, 30,QPen(Qt::darkBlue,2)));
    Lines.append(new Line(  0,-10,  0,-30,QPen(Qt::darkBlue,2)));

    Ports.append(new Port(  0,-30));    
    Ports.append(new Port(-30,  0));
    Ports.append(new Port( 30,  0));
    Ports.append(new Port(  0, 30));
    
    if (elem_type == "cmosinvt")
    {
       Ports.append(new Port( 30,-20));
       Ports.append(new Port( 30, 20));
       Texts.append(new Text( 35, 15,"ref"));
    }
    
    if (elem_type == "opampt")
    {
       Ports.append(new Port( 25, 20));
       Texts.append(new Text( 32, 15,"temp"));
    }

    x1 = -30; y1 = -30;
    x2 =  30; y2 =  30;
    
    tx = x2-8;
    ty = y2;  
  }
  else if ((elem_type == "cmos2nand")   || (elem_type == "cmos2nandx") ||
           (elem_type == "cmos2nandt") || (elem_type == "cmosnand") || (elem_type == "cmosnandt"))
  { 
    Lines.append(new Line(-14,-16,-14, 16,QPen(Qt::darkBlue,2)));
    Lines.append(new Line(-14,-16,  1,-16,QPen(Qt::darkBlue,2)));
    Lines.append(new Line(-14, 16,  1, 16,QPen(Qt::darkBlue,2)));
    
    Arcs.append(new Arc(-12,-16, 26, 32,16*-90,16*180,QPen(Qt::darkBlue,2)));
    Arcs.append(new Arc( 14, -3,  6,  6,  0,16*360,QPen(Qt::darkBlue,2))); 
      
    Lines.append(new Line(-30, 10,-14, 10,QPen(Qt::darkBlue,2)));
    Lines.append(new Line(-30,-10,-14,-10,QPen(Qt::darkBlue,2)));
    Lines.append(new Line( 20,  0, 30,  0,QPen(Qt::darkBlue,2)));
    
    Lines.append(new Line(  0, 16,  0, 30,QPen(Qt::darkBlue,2)));
    Lines.append(new Line(  0,-16,  0,-30,QPen(Qt::darkBlue,2)));

    Ports.append(new Port(  0,-30));
    Ports.append(new Port(-30,-10));    
    Ports.append(new Port(-30, 10));
    Ports.append(new Port( 30,  0));
    Ports.append(new Port(  0, 30));
  
    x1 = -30; y1 = -30;
    x2 =  30; y2 =  30;

    if ((elem_type == "cmos2nandt") || (elem_type == "cmosnandt"))
    {
       Ports.append(new Port( 30,-20));
       Ports.append(new Port( 30, 20));
       Texts.append(new Text( 35, 15,"ref"));
    }
  
    tx = x2-8;
    ty = y2;  
  }
  else if((elem_type == "cmoslatchsrnand") || (elem_type == "cmoslatchsrnandx"))
  {
    Lines.append(new Line(-14,-16,-14, 16,QPen(Qt::darkBlue,2)));
    Lines.append(new Line(-14,-16, 14 ,-16,QPen(Qt::darkBlue,2)));
    Lines.append(new Line(-14, 16, 14, 16,QPen(Qt::darkBlue,2)));
    Lines.append(new Line(14, -16, 14, 16,QPen(Qt::darkBlue,2)));

    //Arcs.append(new Arc( 14, -14,  7,  7,  0,16*360,QPen(Qt::darkBlue,2)));
    Arcs.append(new Arc( 14, 6,  7,  7,  0,16*360,QPen(Qt::darkBlue,2)));
    Arcs.append(new Arc( -21, -14,  7,  7,  0,16*360,QPen(Qt::darkBlue,2)));
    Arcs.append(new Arc( -21, 6,  7,  7,  0,16*360,QPen(Qt::darkBlue,2)));
   
    Lines.append(new Line(-30, 10,-21, 10,QPen(Qt::darkBlue,2)));
    Lines.append(new Line(-30,-10,-21,-10,QPen(Qt::darkBlue,2)));
    Lines.append(new Line(21,10,30,10,QPen(Qt::darkBlue,2)));
    Lines.append(new Line(14,-10,30,-10,QPen(Qt::darkBlue,2)));

    Lines.append(new Line(  0, 16,  0, 30,QPen(Qt::darkBlue,2)));
    Lines.append(new Line(  0,-16,  0,-30,QPen(Qt::darkBlue,2)));

    Ports.append(new Port(  0,-30));
    Ports.append(new Port(-30,-10));
    Ports.append(new Port(-30, 10));
    Ports.append(new Port( 30, 10));
    Ports.append(new Port( 30, -10));
    Ports.append(new Port(  0, 30));

    x1 = -30; y1 = -30;
    x2 =  30; y2 =  30;

   tx = x2-8;
   ty = y2;
  }
  else if ((elem_type == "cmos2nor")   ||
           (elem_type == "cmos2nort") || (elem_type == "cmosnor") || (elem_type == "cmosnort"))
  { 
    //Lines.append(new Line(-14,-16,-14, 16,QPen(Qt::darkBlue,2)));
    Lines.append(new Line(-22,-16,  1,-16,QPen(Qt::darkBlue,2)));
    Lines.append(new Line(-22, 16,  1, 16,QPen(Qt::darkBlue,2)));
    
    Arcs.append(new Arc(-30,-16, 16, 32,16*-90,16*180,QPen(Qt::darkBlue,2)));
    Arcs.append(new Arc( 14, -3,  6,  6,  0,16*360,QPen(Qt::darkBlue,2)));
    
    Arcs.append(new Arc(-18,-16, 32, 32, 0,10*135,QPen(Qt::darkBlue,2))); 
    Arcs.append(new Arc(-18,-16, 32, 32, 0,10*-135,QPen(Qt::darkBlue,2)));
      
    Lines.append(new Line(-30, 10,-15, 10,QPen(Qt::darkBlue,2)));
    Lines.append(new Line(-30,-10,-15,-10,QPen(Qt::darkBlue,2)));
    Lines.append(new Line( 20,  0, 30,  0,QPen(Qt::darkBlue,2)));
    
    Lines.append(new Line(  0, 16,  0, 30,QPen(Qt::darkBlue,2)));
    Lines.append(new Line(  0,-16,  0,-30,QPen(Qt::darkBlue,2)));
    
    Ports.append(new Port(  0,-30));
    Ports.append(new Port(-30,-10));    
    Ports.append(new Port(-30, 10));
    Ports.append(new Port( 30,  0));
    Ports.append(new Port(  0, 30));

    x1 = -30; y1 = -30;
    x2 =  30; y2 =  30;

    if ((elem_type == "cmos2nort") || (elem_type == "cmosnort"))
    {
       Ports.append(new Port( 30,-20));
       Ports.append(new Port( 30, 20));
       Texts.append(new Text( 35, 15,"ref"));
    }
    
    tx = x2-8;
    ty = y2;  
  }
    else if ((elem_type == "mosnbsim3") || (elem_type == "mosnbsim4")  ||
             (elem_type == "mospbsim3") || (elem_type == "mosnbsim3soi5t1") ||
	     (elem_type == "mospbsim3soi5t1"))
  {
     Lines.append(new Line(-14,-13,-14, 13,QPen(Qt::darkBlue,3)));
     Lines.append(new Line(-30,  0,-14,  0,QPen(Qt::darkBlue,2)));

     Lines.append(new Line(-10,-11,  0,-11,QPen(Qt::darkBlue,2)));
     Lines.append(new Line(  0,-11,  0,-30,QPen(Qt::darkBlue,2)));
     Lines.append(new Line(-10, 11,  0, 11,QPen(Qt::darkBlue,2)));
     Lines.append(new Line(  0, 11,  0, 30,QPen(Qt::darkBlue,2)));
     Lines.append(new Line(-10,  0, 20,  0,QPen(Qt::darkBlue,2)));

     Lines.append(new Line(-10,-16,-10, -7,QPen(Qt::darkBlue,3)));
     Lines.append(new Line(-10,  7,-10, 16,QPen(Qt::darkBlue,3)));

     Lines.append(new Line( -4, 24,  4, 20,QPen(Qt::darkBlue,2)));

     // These three lines must be the last.
     Lines.append(new Line(-10, -4,-10,  4,QPen(Qt::darkBlue,3)));
     if ((elem_type == "mospbsim3") || (elem_type == "mospbsim3soi5t1"))
     {
       Lines.append(new Line(  1,  0, -4, -5,QPen(Qt::darkBlue,2)));
       Lines.append(new Line(  1,  0, -4,  5,QPen(Qt::darkBlue,2)));
     }
     else
     {
       Lines.append(new Line( -8,  0, -4, -5,QPen(Qt::darkBlue,2)));
       Lines.append(new Line( -8,  0, -4,  5,QPen(Qt::darkBlue,2)));     
     }

     Ports.append(new Port(  0,-30));
     Ports.append(new Port(-30,  0));
     Ports.append(new Port(  0, 30));
     Ports.append(new Port( 20,  0));
     
     if ((elem_type == "mosnbsim3soi5t1")|| (elem_type == "mospbsim3soi5t1"))
     {
       Texts.append(new Text(20, 15,"Ve"));
       Ports.append(new Port(15, 20));
     }
  
     x1 = -30; y1 = -30;
     x2 =  30; y2 =  30;

     tx = x2+4;
     ty = y1+4;
  }    
  else if ((elem_type == "mosp1") || (elem_type == "mosp2")      ||
           (elem_type == "mosp3") || (elem_type == "mosp9")      ||
	   (elem_type == "mospekv") || (elem_type == "mospekv2")      ||
     (elem_type == "mospekvt"))
  {
     Lines.append(new Line(-14,-13,-14, 13,QPen(Qt::darkBlue,3)));
     Lines.append(new Line(-30,  0,-14,  0,QPen(Qt::darkBlue,2)));

     Lines.append(new Line(-10,-11,  0,-11,QPen(Qt::darkBlue,2)));
     Lines.append(new Line(  0,-11,  0,-30,QPen(Qt::darkBlue,2)));
     Lines.append(new Line(-10, 11,  0, 11,QPen(Qt::darkBlue,2)));
     Lines.append(new Line(  0, 11,  0, 30,QPen(Qt::darkBlue,2)));
     Lines.append(new Line(-10,  0, 20,  0,QPen(Qt::darkBlue,2)));

     Lines.append(new Line(-10,-16,-10, -7,QPen(Qt::darkBlue,3)));
     Lines.append(new Line(-10,  7,-10, 16,QPen(Qt::darkBlue,3)));

     Lines.append(new Line( -4, 24,  4, 20,QPen(Qt::darkBlue,2)));

     // These three lines must be the last.
     Lines.append(new Line(-10, -4,-10,  4,QPen(Qt::darkBlue,3)));
     Lines.append(new Line(  1,  0, -4, -5,QPen(Qt::darkBlue,2)));
     Lines.append(new Line(  1,  0, -4,  5,QPen(Qt::darkBlue,2)));

     Ports.append(new Port(-30,  0));
     Ports.append(new Port(  0,-30));
     Ports.append(new Port(  0, 30));
     Ports.append(new Port( 20,  0));

     if (elem_type == "mospekvt")
     {
      Ports.append(new Port( 10,-20));
      Ports.append(new Port( 10, 20));
      Texts.append(new Text( 15, 15,"ref"));
     }
       
     x1 = -30; y1 = -30;
     x2 =  30; y2 =  30;

     tx = x2+4;
     ty = y1+4;
  }    
  else if ((elem_type == "mosn1") || (elem_type == "mosn2")         || 
           (elem_type == "mosn3") || (elem_type == "mosn9")         ||
	   (elem_type == "mosnekv")   || (elem_type == "mosnekv2")  || 
	   (elem_type == "mosnekv3")  || (elem_type == "mosnldmet") ||
	   (elem_type == "mosntftq") || (elem_type == "mosnekvt"))
  {
     Lines.append(new Line(-14,-13,-14, 13,QPen(Qt::darkBlue,3)));
     Lines.append(new Line(-30,  0,-14,  0,QPen(Qt::darkBlue,2)));

     Lines.append(new Line(-10,-11,  0,-11,QPen(Qt::darkBlue,2)));
     Lines.append(new Line(  0,-11,  0,-30,QPen(Qt::darkBlue,2)));
     Lines.append(new Line(-10, 11,  0, 11,QPen(Qt::darkBlue,2)));
     Lines.append(new Line(  0, 11,  0, 30,QPen(Qt::darkBlue,2)));
     Lines.append(new Line(-10,  0, 20,  0,QPen(Qt::darkBlue,2)));

     Lines.append(new Line(-10,-16,-10, -7,QPen(Qt::darkBlue,3)));
     Lines.append(new Line(-10,  7,-10, 16,QPen(Qt::darkBlue,3)));

     Lines.append(new Line( -4, 24,  4, 20,QPen(Qt::darkBlue,2)));

     // These three lines must be the last.
     Lines.append(new Line(-10, -4,-10,  4,QPen(Qt::darkBlue,3)));
     Lines.append(new Line( -9,  0, -4, -5,QPen(Qt::darkBlue,2)));
     Lines.append(new Line( -9,  0, -4,  5,QPen(Qt::darkBlue,2)));

     Ports.append(new Port(-30,  0));
     Ports.append(new Port(  0,-30));
     Ports.append(new Port(  0, 30));
     Ports.append(new Port( 20,  0));
     
     if ((elem_type == "mosnldmet")|| (elem_type == "mosnekvt"))
     {
       Ports.append(new Port( 10,-20));
       Ports.append(new Port( 10, 20));
       Texts.append(new Text( 15, 15,"ref"));    
     }
       
     x1 = -30; y1 = -30;
     x2 =  30; y2 =  30;

     tx = x2+4;
     ty = y1+4;
  }
  else if ((elem_type == "nport") && (nport_str == "nport4"))
  {
    Lines.append(new Line(-14,-14, 14,-14,QPen(Qt::darkBlue,2)));
    Lines.append(new Line(-14, 14, 14, 14,QPen(Qt::darkBlue,2)));
    Lines.append(new Line(-14,-14,-14, 14,QPen(Qt::darkBlue,2)));
    Lines.append(new Line( 14,-14, 14, 14,QPen(Qt::darkBlue,2)));

    Lines.append(new Line(-30, 10,-14, 10,QPen(Qt::darkBlue,2)));
    Lines.append(new Line(-30,-10,-14,-10,QPen(Qt::darkBlue,2)));    
    Lines.append(new Line( 14,-10, 30,-10,QPen(Qt::darkBlue,2)));
    Lines.append(new Line( 14, 10, 30, 10,QPen(Qt::darkBlue,2)));

    Ports.append(new Port(-30,-10));
    Ports.append(new Port( 30,-10));
    Ports.append(new Port(-30, 10));
    Ports.append(new Port( 30, 10));
    
    x1 = -30; y1 = -17;
    x2 =  30; y2 =  17;

    tx = x1+4;
    ty = y2+4;
    
    elem_type = "nport4"; 
  } 
  else if ((elem_type == "nport") && (nport_str == "nport6"))
  {
    Lines.append(new Line(-14,-14, 14,-14,QPen(Qt::darkBlue,2)));
    Lines.append(new Line(-14, 14, 14, 14,QPen(Qt::darkBlue,2)));
    Lines.append(new Line(-14,-14,-14, 14,QPen(Qt::darkBlue,2)));
    Lines.append(new Line( 14,-14, 14, 14,QPen(Qt::darkBlue,2)));

    Lines.append(new Line(-30, 10,-14, 10,QPen(Qt::darkBlue,2)));
    Lines.append(new Line(-30,-10,-14,-10,QPen(Qt::darkBlue,2)));    
    Lines.append(new Line(-30,  0,-14,  0,QPen(Qt::darkBlue,2)));    
    Lines.append(new Line( 14,  0, 30,  0,QPen(Qt::darkBlue,2)));
    Lines.append(new Line( 14,-10, 30,-10,QPen(Qt::darkBlue,2)));
    Lines.append(new Line( 14, 10, 30, 10,QPen(Qt::darkBlue,2)));

    Ports.append(new Port(-30,-10));
    Ports.append(new Port( 30,-10));
    Ports.append(new Port(-30,  0));
    Ports.append(new Port( 30,  0));
    Ports.append(new Port(-30, 10));
    Ports.append(new Port( 30, 10));
    
    x1 = -30; y1 = -17;
    x2 =  30; y2 =  17;

    tx = x1+4;
    ty = y2+4;
    
    elem_type = "nport6";
  }     
  else if ((elem_type == "nport") && (nport_str == "nport8"))        
  {

    Lines.append(new Line(-14,-24, 14,-24,QPen(Qt::darkBlue,2)));
    Lines.append(new Line(-14, 14, 14, 14,QPen(Qt::darkBlue,2)));
    Lines.append(new Line(-14,-24,-14, 14,QPen(Qt::darkBlue,2)));
    Lines.append(new Line( 14,-24, 14, 14,QPen(Qt::darkBlue,2)));

    Lines.append(new Line(-30,-20,-14,-20,QPen(Qt::darkBlue,2)));
    Lines.append(new Line(-30,-10,-14,-10,QPen(Qt::darkBlue,2)));    
    Lines.append(new Line(-30,  0,-14,  0,QPen(Qt::darkBlue,2)));
    Lines.append(new Line(-30, 10,-14, 10,QPen(Qt::darkBlue,2)));    

    Lines.append(new Line( 14,-20, 30,-20,QPen(Qt::darkBlue,2)));
    Lines.append(new Line( 14,-10, 30,-10,QPen(Qt::darkBlue,2)));
    Lines.append(new Line( 14,  0, 30,  0,QPen(Qt::darkBlue,2)));
    Lines.append(new Line( 14, 10, 30, 10,QPen(Qt::darkBlue,2)));

    Ports.append(new Port(-30,-20));
    Ports.append(new Port( 30,-20));
    
    Ports.append(new Port(-30,-10));
    Ports.append(new Port( 30,-10));

    Ports.append(new Port(-30,  0));
    Ports.append(new Port( 30,  0));

    Ports.append(new Port(-30, 10));
    Ports.append(new Port( 30, 10));
    
    x1 = -30; y1 = -27;
    x2 =  30; y2 =  17;

    tx = x1+4;
    ty = y2+4;
    
    elem_type = "nport8";  
  }
  else if ((elem_type == "nport") && (nport_str == "nport10"))        
  {
    Lines.append(new Line(-14,-24, 14,-24,QPen(Qt::darkBlue,2)));
    Lines.append(new Line(-14, 24, 14, 24,QPen(Qt::darkBlue,2)));
    Lines.append(new Line(-14,-24,-14, 24,QPen(Qt::darkBlue,2)));
    Lines.append(new Line( 14,-24, 14, 24,QPen(Qt::darkBlue,2)));

    Lines.append(new Line(-30,-20,-14,-20,QPen(Qt::darkBlue,2)));
    Lines.append(new Line(-30,-10,-14,-10,QPen(Qt::darkBlue,2)));    
    Lines.append(new Line(-30,  0,-14,  0,QPen(Qt::darkBlue,2)));
    Lines.append(new Line(-30, 10,-14, 10,QPen(Qt::darkBlue,2)));    
    Lines.append(new Line(-30, 20,-14, 20,QPen(Qt::darkBlue,2)));    

    Lines.append(new Line( 14,-20, 30,-20,QPen(Qt::darkBlue,2)));
    Lines.append(new Line( 14,-10, 30,-10,QPen(Qt::darkBlue,2)));
    Lines.append(new Line( 14,  0, 30,  0,QPen(Qt::darkBlue,2)));
    Lines.append(new Line( 14, 10, 30, 10,QPen(Qt::darkBlue,2)));
    Lines.append(new Line( 14, 20, 30, 20,QPen(Qt::darkBlue,2)));

    Ports.append(new Port(-30,-20));
    Ports.append(new Port( 30,-20));
    
    Ports.append(new Port(-30,-10));
    Ports.append(new Port( 30,-10));

    Ports.append(new Port(-30,  0));
    Ports.append(new Port( 30,  0));

    Ports.append(new Port(-30, 10));
    Ports.append(new Port( 30, 10));

    Ports.append(new Port(-30, 20));
    Ports.append(new Port( 30, 20));
    
    x1 = -30; y1 = -27;
    x2 =  30; y2 =  27;

    tx = x1+4;
    ty = y2+4;
    
    elem_type = "nport10"; 
  }    
  else if (elem_type == "open")
  {
    Lines.append(new Line(-30,  0,-12,  0,QPen(Qt::darkBlue,2)));
    Lines.append(new Line( 12,  0, 30,  0,QPen(Qt::darkBlue,2)));
    
    Ports.append(new Port(-30,  0));
    Ports.append(new Port( 30,  0));

    x1 = -30; y1 = -11;
    x2 =  30; y2 =  11;

    tx = x1+10;
    ty = y2-11;
 
  }
  else if ((elem_type == "r")               ||
           (elem_type == "res")             ||
	   (elem_type == "resistor")        ||
	   (elem_type == "resistorphyn")    ||
	   (elem_type == "resistorphyp")    ||
	   (elem_type == "resistorphypoly") ||
	   (elem_type == "resistort")) 
  {
    Lines.append(new Line(-30,  0,-18,  0,QPen(Qt::darkBlue,2)));
    Lines.append(new Line(-18,  0,-15, -7,QPen(Qt::darkBlue,2)));
    Lines.append(new Line(-15, -7, -9,  7,QPen(Qt::darkBlue,2)));
    Lines.append(new Line( -9,  7, -3, -7,QPen(Qt::darkBlue,2)));
    Lines.append(new Line( -3, -7,  3,  7,QPen(Qt::darkBlue,2)));
    Lines.append(new Line(  3,  7,  9, -7,QPen(Qt::darkBlue,2)));
    Lines.append(new Line(  9, -7, 15,  7,QPen(Qt::darkBlue,2)));
    Lines.append(new Line( 15,  7, 18,  0,QPen(Qt::darkBlue,2)));
    Lines.append(new Line( 18,  0, 30,  0,QPen(Qt::darkBlue,2)));
    
    Ports.append(new Port(-30,  0));
    Ports.append(new Port( 30,  0));

    x1 = -30; y1 = -11;
    x2 =  30; y2 =  11;

    if ((elem_type == "resistorphyn") ||	   
        (elem_type == "resistorphyp") ||
        (elem_type == "resistorphypoly")) 

    {
       Lines.append(new Line(  0,  0,  0, 20,QPen(Qt::darkBlue,2)));
       Ports.append(new Port(  0, 20));

       x1 = -30; y1 = -11;
       x2 =  30; y2 =  20;     
    }

    if (elem_type == "resistort")  

    {
       Ports.append(new Port( 15, 15));
       Texts.append(new Text( 17, 14,"Tout"));

       Ports.append(new Port(-15, 15));
       Texts.append(new Text( -7, 14,"Tref"));

       x1 = -30; y1 = -11;
       x2 =  30; y2 =  22;     
    }

    tx = x1+4;
    ty = y2+4;
  }
  
  else if ((elem_type == "thermalblockbc1")        ||
           (elem_type == "thermalheatsink1")       ||
           (elem_type == "thermalheatsinkmmic1")   ||
	        (elem_type == "thermaltransf")           ||
          (elem_type == "thermaltest"))
  {

    Lines.append(new Line(-14,-14, 14,-14,QPen(Qt::darkBlue,2)));
    Lines.append(new Line(-14, 14, 14, 14,QPen(Qt::darkBlue,2)));
    Lines.append(new Line(-14,-14,-14, 14,QPen(Qt::darkBlue,2)));
    Lines.append(new Line( 14,-14, 14, 14,QPen(Qt::darkBlue,2)));

    Lines.append(new Line(-30,  0,-14,  0,QPen(Qt::darkBlue,2)));
    Lines.append(new Line( 30,  0, 14,  0,QPen(Qt::darkBlue,2)));    

    Ports.append(new Port(-30,  0));
    Ports.append(new Port( 30,  0));    
 
    x1 = -30; y1 = -17;
    x2 =  30; y2 =  17;

    tx = x1+4;
    ty = y2+4;
      
  }   

  else if (elem_type == "thermalint")
  {

    Lines.append(new Line(-14,-14, 14,-14,QPen(Qt::darkBlue,2)));
    Lines.append(new Line(-14, 14, 14, 14,QPen(Qt::darkBlue,2)));
    Lines.append(new Line(-14,-14,-14, 14,QPen(Qt::darkBlue,2)));
    Lines.append(new Line( 14,-14, 14, 14,QPen(Qt::darkBlue,2)));

    Lines.append(new Line(-30,  0,-14,  0,QPen(Qt::darkBlue,2)));
    Lines.append(new Line( 30,  0, 14,  0,QPen(Qt::darkBlue,2)));    
    Lines.append(new Line(  0, 30,  0, 14,QPen(Qt::darkBlue,2)));

    Ports.append(new Port(-30,  0));
    Ports.append(new Port( 30,  0));    
    Ports.append(new Port(  0, 30));    
 
    x1 = -30; y1 = -17;
    x2 =  30; y2 =  30;

    tx = x2+4;
    ty = y2+4;
      
  }   
  else if (elem_type == "thermalshunt")
  {

    Lines.append(new Line(-14,-14,  5,-14,QPen(Qt::darkBlue,2)));
    Lines.append(new Line(-14, 14,  5, 14,QPen(Qt::darkBlue,2)));
    Lines.append(new Line(-14,-14,-14, 14,QPen(Qt::darkBlue,2)));
    Lines.append(new Line(  5,-14,  5, 14,QPen(Qt::darkBlue,2)));

    Lines.append(new Line(-30, 10,-14, 10,QPen(Qt::darkBlue,2)));
    Lines.append(new Line(-30,-10,-14,-10,QPen(Qt::darkBlue,2)));    
    Lines.append(new Line(-30,  0,-14,  0,QPen(Qt::darkBlue,2)));

    Ports.append(new Port(-30,-10));
    Ports.append(new Port(-30,  0));    
    Ports.append(new Port(-30, 10));    
 
    x1 = -30; y1 = -17;
    x2 =   7; y2 =  17;

    tx = x1+4;
    ty = y2+4;
  
  }   
  else if (elem_type == "tlinp4")
  {
    Lines.append(new Line(-30,-12,-16,-12,QPen(Qt::darkBlue,2)));
    Lines.append(new Line(-30,-30,-30,-12,QPen(Qt::darkBlue,2)));
    Lines.append(new Line( 20,-12, 30,-12,QPen(Qt::darkBlue,2)));
    Lines.append(new Line( 30,-30, 30,-12,QPen(Qt::darkBlue,2)));
    Lines.append(new Line(-11,-20, 25,-20,QPen(Qt::darkBlue,2)));
    Lines.append(new Line(-21, -4, 15, -4,QPen(Qt::darkBlue,2)));
    Lines.append(new Line(-11,-20,-21, -4,QPen(Qt::darkBlue,2)));
    Lines.append(new Line( 25,-20, 15, -4,QPen(Qt::darkBlue,2)));

    Lines.append(new Line(-30, 12,-20, 12,QPen(Qt::darkBlue,2)));
    Lines.append(new Line(-30, 30,-30, 12,QPen(Qt::darkBlue,2)));
    Lines.append(new Line( 16, 12, 30, 12,QPen(Qt::darkBlue,2)));
    Lines.append(new Line( 30, 30, 30, 12,QPen(Qt::darkBlue,2)));
    Lines.append(new Line(-15,  4, 21,  4,QPen(Qt::darkBlue,2)));
    Lines.append(new Line(-25, 20, 11, 20,QPen(Qt::darkBlue,2)));
    Lines.append(new Line(-15,  4,-25, 20,QPen(Qt::darkBlue,2)));
    Lines.append(new Line( 21,  4, 11, 20,QPen(Qt::darkBlue,2)));

    Lines.append(new Line(-22,-16,-26, -8,QPen(Qt::darkBlue,2)));

    Ports.append(new Port(-30,-30));
    Ports.append(new Port(-30, 30));
    Ports.append(new Port( 30,-30));
    Ports.append(new Port( 30, 30));


    x1 = -30; y1 =-33;
    x2 =  30; y2 = 33;

    tx = x1+4;
    ty = y2+4;

  }
  else if (elem_type == "vcvs")
  {
    Arcs.append(new Arc(0,-11, 23, 23,  0, 16*360,QPen(Qt::darkBlue,2)));

    Lines.append(new Line(-30,-30,-12,-30,QPen(Qt::darkBlue,2)));
    Lines.append(new Line(-30, 30,-12, 30,QPen(Qt::darkBlue,2)));
    Lines.append(new Line( 11,-30, 30,-30,QPen(Qt::darkBlue,2)));
    Lines.append(new Line( 11, 30, 30, 30,QPen(Qt::darkBlue,2)));

    Lines.append(new Line(-12,-30,-12,-23,QPen(Qt::darkBlue,2)));
    Lines.append(new Line(-12, 30,-12, 23,QPen(Qt::darkBlue,2)));
    Lines.append(new Line( 11,-30, 11,-11,QPen(Qt::darkBlue,2)));
    Lines.append(new Line( 11, 30, 11, 11,QPen(Qt::darkBlue,2)));

    Lines.append(new Line(-12,-18,-12, 18,QPen(Qt::darkBlue,2)));
    Lines.append(new Line(-12, 18,-17,  9,QPen(Qt::darkBlue,2)));
    Lines.append(new Line(-12, 18, -8,  9,QPen(Qt::darkBlue,2)));

    Lines.append(new Line( 19,-21, 19,-15,QPen(Qt::red,1)));
    Lines.append(new Line( 16,-18, 22,-18,QPen(Qt::red,1)));
    Lines.append(new Line( 16, 18, 22, 18,QPen(Qt::black,1)));

    Lines.append(new Line(-25,-27, 25,-27,QPen(Qt::darkGray,1)));
    Lines.append(new Line( 25,-27, 25, 27,QPen(Qt::darkGray,1)));
    Lines.append(new Line( 25, 27,-25, 27,QPen(Qt::darkGray,1)));
    Lines.append(new Line(-25, 27,-25,-27,QPen(Qt::darkGray,1)));

    Ports.append(new Port( 30,-30));
    Ports.append(new Port( 30, 30));
    Ports.append(new Port(-30,-30));
    Ports.append(new Port(-30, 30));


    x1 = -30; y1 = -30;
    x2 =  30; y2 =  30;

    tx = x1+4;
    ty = y2+4;  
  }

  else if ((elem_type == "vvtb_cann") || (elem_type == "vvtb_tanh"))
  {
    Arcs.append(new Arc(0,-11, 23, 23,  0, 16*360,QPen(Qt::darkBlue,2)));

    Lines.append(new Line(-30,-30,-12,-30,QPen(Qt::darkBlue,2)));
    Lines.append(new Line(-30, 30,-12, 30,QPen(Qt::darkBlue,2)));
    Lines.append(new Line( 11,-30, 30,-30,QPen(Qt::darkBlue,2)));
    Lines.append(new Line( 11, 30, 30, 30,QPen(Qt::darkBlue,2)));

    Lines.append(new Line(-12,-30,-12,-23,QPen(Qt::darkBlue,2)));
    Lines.append(new Line(-12, 30,-12, 23,QPen(Qt::darkBlue,2)));
    Lines.append(new Line( 11,-30, 11,-11,QPen(Qt::darkBlue,2)));
    Lines.append(new Line( 11, 30, 11, 11,QPen(Qt::darkBlue,2)));

    Arcs.append(new Arc(-16,-23, 8, 8,  0, 16*360,QPen(Qt::darkBlue,2)));
    Arcs.append(new Arc(-16, 15, 8, 8,  0, 16*360,QPen(Qt::darkBlue,2)));


    Lines.append(new Line( 19,-21, 19,-15,QPen(Qt::red,1)));
    Lines.append(new Line( 16,-18, 22,-18,QPen(Qt::red,1)));
    Lines.append(new Line( 16, 18, 22, 18,QPen(Qt::black,1)));

    Lines.append(new Line(-25,-27, 25,-27,QPen(Qt::darkGray,1)));
    Lines.append(new Line( 25,-27, 25, 27,QPen(Qt::darkGray,1)));
    Lines.append(new Line( 25, 27,-25, 27,QPen(Qt::darkGray,1)));
    Lines.append(new Line(-25, 27,-25,-27,QPen(Qt::darkGray,1)));

    Ports.append(new Port( 30,-30));
    Ports.append(new Port( 30, 30));
    Ports.append(new Port(-30,-30));
    Ports.append(new Port(-30, 30));

    x1 = -30; y1 = -30;
    x2 =  30; y2 =  30;

    tx = x1+4;
    ty = y2+4;  
    }


  else if ((elem_type == "vccs") || (elem_type == "vccspoly") || (elem_type == "cccs"))
  {
    Arcs.append(new Arc(0,-11, 23, 23,  0, 16*360,QPen(Qt::darkBlue,2)));
    Lines.append(new Line( 11, -7, 11,  7,QPen(Qt::darkBlue,3)));
    Lines.append(new Line( 11,  6, 15,  1,QPen(Qt::darkBlue,3)));
    Lines.append(new Line( 11,  6,  7,  1,QPen(Qt::darkBlue,3)));

    Lines.append(new Line(-30,-30,-12,-30,QPen(Qt::darkBlue,2)));
    Lines.append(new Line(-30, 30,-12, 30,QPen(Qt::darkBlue,2)));
    Lines.append(new Line( 11,-30, 30,-30,QPen(Qt::darkBlue,2)));
    Lines.append(new Line( 11, 30, 30, 30,QPen(Qt::darkBlue,2)));

    Lines.append(new Line(-12,-30,-12,-23,QPen(Qt::darkBlue,2)));
    Lines.append(new Line(-12, 30,-12, 23,QPen(Qt::darkBlue,2)));
    Lines.append(new Line( 11,-30, 11,-11,QPen(Qt::darkBlue,2)));
    Lines.append(new Line( 11, 30, 11, 11,QPen(Qt::darkBlue,2)));

    Lines.append(new Line(-12,-18,-12, 18,QPen(Qt::darkBlue,2)));
    Lines.append(new Line(-12, 18,-17,  9,QPen(Qt::darkBlue,2)));
    Lines.append(new Line(-12, 18, -8,  9,QPen(Qt::darkBlue,2)));

    Lines.append(new Line(-25,-27, 25,-27,QPen(Qt::darkGray,1)));
    Lines.append(new Line( 25,-27, 25, 27,QPen(Qt::darkGray,1)));
    Lines.append(new Line( 25, 27,-25, 27,QPen(Qt::darkGray,1)));
    Lines.append(new Line(-25, 27,-25,-27,QPen(Qt::darkGray,1)));

    Ports.append(new Port( 30,-30));
    Ports.append(new Port( 30, 30));
    Ports.append(new Port(-30,-30));
    Ports.append(new Port(-30, 30));

    x1 = -30; y1 = -30;
    x2 =  30; y2 =  30;

    tx = x1+4;
    ty = y2+4;
  }
  else if (elem_type == "vccsd")
  {
    Arcs.append(new Arc(0,-11, 23, 23,  0, 16*360,QPen(Qt::darkBlue,2)));
    Lines.append(new Line( 11, -7, 11,  7,QPen(Qt::darkBlue,3)));
    Lines.append(new Line( 11,  6, 15,  1,QPen(Qt::darkBlue,3)));
    Lines.append(new Line( 11,  6,  7,  1,QPen(Qt::darkBlue,3)));

    Lines.append(new Line(-30,-30,-12,-30,QPen(Qt::darkBlue,2)));
    Lines.append(new Line(-30, 30,-12, 30,QPen(Qt::darkBlue,2)));
    Lines.append(new Line( 11,-30, 30,-30,QPen(Qt::darkBlue,2)));
    Lines.append(new Line( 11, 30, 30, 30,QPen(Qt::darkBlue,2)));

    Lines.append(new Line(-12,-30,-12,-23,QPen(Qt::darkBlue,2)));
    Lines.append(new Line(-12, 30,-12, 23,QPen(Qt::darkBlue,2)));
    Lines.append(new Line( 11,-30, 11,-11,QPen(Qt::darkBlue,2)));
    Lines.append(new Line( 11, 30, 11, 11,QPen(Qt::darkBlue,2)));
    
    Arcs.append(new Arc(-16,-23, 8, 8,  0, 16*360,QPen(Qt::darkBlue,2)));
    Arcs.append(new Arc(-16, 15, 8, 8,  0, 16*360,QPen(Qt::darkBlue,2)));

    Lines.append(new Line(-25,-27, 25,-27,QPen(Qt::darkGray,1)));
    Lines.append(new Line( 25,-27, 25, 27,QPen(Qt::darkGray,1)));
    Lines.append(new Line( 25, 27,-25, 27,QPen(Qt::darkGray,1)));
    Lines.append(new Line(-25, 27,-25,-27,QPen(Qt::darkGray,1)));

    Ports.append(new Port( 30,-30));
    Ports.append(new Port( 30, 30));
    Ports.append(new Port(-30,-30));
    Ports.append(new Port(-30, 30));

      x1 = -30; y1 = -30;
      x2 =  30; y2 =  30;

      tx = x1+4;
      ty = y2+4;
  }

  else if ((elem_type == "vct")|| (elem_type == "vctb_cann") || (elem_type == "vctb_tanh"))
  {
    Arcs.append(new Arc(0,-11, 23, 23,  0, 16*360,QPen(Qt::darkBlue,2)));
    Lines.append(new Line( 11, -7, 11,  7,QPen(Qt::darkBlue,3)));
    Lines.append(new Line( 11,  6, 15,  1,QPen(Qt::darkBlue,3)));
    Lines.append(new Line( 11,  6,  7,  1,QPen(Qt::darkBlue,3)));

    Lines.append(new Line(-30,-30,-12,-30,QPen(Qt::darkBlue,2)));
    Lines.append(new Line(-30, 30,-12, 30,QPen(Qt::darkBlue,2)));
    Lines.append(new Line( 11,-30, 30,-30,QPen(Qt::darkBlue,2)));
    Lines.append(new Line( 11, 30, 30, 30,QPen(Qt::darkBlue,2)));

    Lines.append(new Line(-12,-30,-12,-23,QPen(Qt::darkBlue,2)));
    Lines.append(new Line(-12, 30,-12, 23,QPen(Qt::darkBlue,2)));
    Lines.append(new Line( 11,-30, 11,-11,QPen(Qt::darkBlue,2)));
    Lines.append(new Line( 11, 30, 11, 11,QPen(Qt::darkBlue,2)));

    Arcs.append(new Arc(-16,-23, 8, 8,  0, 16*360,QPen(Qt::darkBlue,2)));
    Arcs.append(new Arc(-16, 15, 8, 8,  0, 16*360,QPen(Qt::darkBlue,2)));
    
    Lines.append(new Line(-25,-27, 25,-27,QPen(Qt::darkGray,1)));
    Lines.append(new Line( 25,-27, 25, 27,QPen(Qt::darkGray,1)));
    Lines.append(new Line( 25, 27,-25, 27,QPen(Qt::darkGray,1)));
    Lines.append(new Line(-25, 27,-25,-27,QPen(Qt::darkGray,1)));

    Ports.append(new Port(-30,-30));
    Ports.append(new Port(-30, 30));    
    Ports.append(new Port( 30,-30));
    Ports.append(new Port( 30, 30));


    x1 = -30; y1 = -30;
    x2 =  30; y2 =  30;

    tx = x1+4;
    ty = y2+4;
  
  }
  else if (elem_type == "vexp")
  {
    Arcs.append(new Arc(-12,-12, 25, 25, 0, 16*360,QPen(Qt::darkBlue,2)));
    Lines.append(new Line(-30,  0,-12,  0,QPen(Qt::darkBlue,2)));
    Lines.append(new Line( 30,  0, 12,  0,QPen(Qt::darkBlue,2)));
    Lines.append(new Line( 18,  5, 18, 11,QPen(Qt::red,1)));
    Lines.append(new Line( 21,  8, 15,  8,QPen(Qt::red,1)));
    Lines.append(new Line(-18,  5,-18, 11,QPen(Qt::black,1)));

    Arcs.append(new Arc( -12, -6, 20, 35,16*10,16*100,QPen(Qt::darkBlue,2)));

    Ports.append(new Port( 30,  0));
    Ports.append(new Port(-30,  0));    


    x1 = -30; y1 = -14;
    x2 =  30; y2 =  14;

    tx = x1+4;
    ty = y2+4;
  
  }
  else if ((elem_type == "vpulse") || (elem_type == "vpulsexp") || (elem_type == "vpulseag") || (elem_type == "vpulserf") || (elem_type == "vlfmpulse"))
  {
    Arcs.append(new Arc(-12,-12, 25, 25,     0, 16*360,QPen(Qt::darkBlue,2)));
    Lines.append(new Line(-30,  0,-12,  0,QPen(Qt::darkBlue,2)));
    Lines.append(new Line( 30,  0, 12,  0,QPen(Qt::darkBlue,2)));
    Lines.append(new Line( 18,  5, 18, 11,QPen(Qt::red,1)));
    Lines.append(new Line( 21,  8, 15,  8,QPen(Qt::red,1)));
    Lines.append(new Line(-18,  5,-18, 11,QPen(Qt::black,1)));

    Lines.append(new Line(  6, -3,  6,  3,QPen(Qt::darkBlue,2)));
    Lines.append(new Line( -6, -7, -6, -3,QPen(Qt::darkBlue,2)));
    Lines.append(new Line( -6,  3, -6,  7,QPen(Qt::darkBlue,2)));
    Lines.append(new Line( -6, -3,  6, -3,QPen(Qt::darkBlue,2)));
    Lines.append(new Line( -6,  3,  6,  3,QPen(Qt::darkBlue,2)));

    Ports.append(new Port( 30,  0));
    Ports.append(new Port(-30,  0));


    x1 = -30; y1 = -14;
    x2 =  30; y2 =  14;

    tx = x1+4;
    ty = y2+4;
  }
  else if ((elem_type == "vsource") || (elem_type == "vam") ||
           (elem_type == "vsffm") || (elem_type == "vtwotone") || (elem_type == "vwnoise"))
  {
    Arcs.append(new Arc(-12,-12, 25, 25,     0, 16*360,QPen(Qt::darkBlue,2)));
    Arcs.append(new Arc( -3, -7,  8,  8,16*270, 16*180,QPen(Qt::darkBlue,2)));
    Arcs.append(new Arc( -3,  0,  8,  8, 16*90, 16*180,QPen(Qt::darkBlue,2)));
    Lines.append(new Line(-30,  0,-12,  0,QPen(Qt::darkBlue,2)));
    Lines.append(new Line( 30,  0, 12,  0,QPen(Qt::darkBlue,2)));
    Lines.append(new Line( 18,  5, 18, 11,QPen(Qt::red,1)));
    Lines.append(new Line( 21,  8, 15,  8,QPen(Qt::red,1)));
    Lines.append(new Line(-18,  5,-18, 11,QPen(Qt::black,1)));

    Ports.append(new Port( 30,  0));
    Ports.append(new Port(-30,  0));


    x1 = -30; y1 = -14;
    x2 =  30; y2 =  14;

    tx = x1+4;
    ty = y2+4;
  }
  else
  {//YOUR ELEMENT WAS NOT FOUND SO YOU NEED TO ADD IT ABOVE
     QMessageBox::critical(0, QObject::tr("Error"),
         QObject::tr("The code to draw your element doesnt exist. Reference interface.cpp"));
  }
}

Component* Interface::newOne()
{
  //QString element_type = elem_type.c_str();
  // Added by Shivam
   QString element_type = elem_type;
  return new Interface(element_type);
}

Component* Interface::info(QString& Name, char* &BitmapFile, bool getNewOne)
{
  //const char* name = Name.data();
  // Changed by Shivam
  const QChar* name = Name.data();
  BitmapFile = (char*)name;
  if(getNewOne)  return new Interface(Name);
  return 0;
}

