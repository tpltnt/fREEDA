/***************************************************************************
                          plot.cpp  -  description
                             -------------------
 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#include "plot.h"
#include "guimain.h"


Plot::Plot()
{
  Description = QObject::tr("Plot Terminal");

  QFont f = QucsSettings.font;
  f.setWeight(QFont::Light);
  f.setPointSizeFloat(12.0);
  QFontMetrics  metrics(f);
  QSize r = metrics.size(0, QObject::tr("Plot Terminal"));
  int xb = r.width()  >> 1;
  int yb = r.height() >> 1;

  Lines.append(new Line(-xb, -yb, -xb,  yb,QPen(Qt::darkBlue,2)));
  Lines.append(new Line(-xb,  yb,  xb+3,yb,QPen(Qt::darkBlue,2)));
  Texts.append(new Text(-xb+4,  -yb-3, QObject::tr("Plot Terminal"),
			QColor(0,0,0), 12.0));

  x1 = -xb-3;  y1 = -yb-5;
  x2 =  xb+9; y2 =  yb+3;

  tx = x1+4;
  ty = y2+4;
  Model = ".out";
  Name  = "plot";

  Props.append(new Property("Mode", "plot", false));
  Props.append(new Property("Qualifier", "term", true,
  			QObject::tr("Qualifier type.(term, element)")));
  Props.append(new Property("Terminal", "", true,
  			QObject::tr("Name of terminal collecting data")));
  Props.append(new Property("Current Element", "", false,
  			QObject::tr("Current element for plotting.")));
  Props.append(new Property("Current Terminal", "", false,
  			QObject::tr("Name of current terminal collecting data")));			
  Props.append(new Property("Network Operator", "vf", true,
  			QObject::tr("Network Operations(vf, if, xf, vt, it, ut)")));			
  Props.append(new Property("Qualifier", "", false,
  			QObject::tr("Qualifier type.(term, element)")));			
  Props.append(new Property("Terminal", "", false,
  			QObject::tr("Name of terminal collecting data.")));	
  Props.append(new Property("Current Element", "", false,
  			QObject::tr("Current element for plotting.")));
  Props.append(new Property("Current Terminal", "", false,
  			QObject::tr("Name of current terminal collecting data")));					
  Props.append(new Property("Network Operator", "", false,
  			QObject::tr("Network Operations(vf, if, xf, vt, it, ut)")));			
  Props.append(new Property("Operators", "", false,
  			QObject::tr("add, sub.. db, db10.. smpltime, sweepfrq.. ")));
  Props.append(new Property("Output file", "test.out", true,
  			QObject::tr("Name of output file.")));

}

Plot::~Plot()
{
}

Component* Plot::newOne()
{
  return new Plot();
}

Component* Plot::info(QString& Name, char* &BitmapFile, bool getNewOne)
{
  Name = QObject::tr("Plot Terminal");
  BitmapFile = "plot";

  if(getNewOne)  return new Plot();
  return 0;
}
