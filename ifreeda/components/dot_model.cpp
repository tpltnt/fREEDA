/***************************************************************************
                          model.cpp  -  description
***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#include "dot_model.h"
#include "guimain.h"


Dot_Model::Dot_Model()
{
  Description = QObject::tr("Model File");

  QFont f = QucsSettings.font;
  f.setWeight(QFont::Light);
  f.setPointSizeFloat(12.0);
  QFontMetrics  metrics(f);
  QSize r = metrics.size(0, QObject::tr("Model File"));
  int xb = r.width()  >> 1;
  int yb = r.height() >> 1;

  Lines.append(new Line(-xb, -yb, -xb,  yb,QPen(Qt::darkBlue,2)));
  Lines.append(new Line(-xb,  yb,  xb+3,yb,QPen(Qt::darkBlue,2)));
  Texts.append(new Text(-xb+4,  -yb-3, QObject::tr("Model File"),
			QColor(0,0,0), 12.0));

  x1 = -xb-3;  y1 = -yb-5;
  x2 =  xb+9; y2 =  yb+3;

  tx = x1+4;
  ty = y2+4;
  Model = ".model";
  Name  = "model";

/*  Props.append(new Property("Model name", "", true,
  			QObject::tr("Name of the model.")));
  Props.append(new Property("Element name", "", true,
  			QObject::tr("Name of the element that will use the model.")));  */
			
			
  Props.append(new Property("Model file", "", true,
  			QObject::tr("Name of model file that contains all the model parameters.")));
}

Dot_Model::~Dot_Model()
{
}

Component* Dot_Model::newOne()
{
  return new Dot_Model();
}

Component* Dot_Model::info(QString& Name, char* &BitmapFile, bool getNewOne)
{
  Name = QObject::tr("Model File");
  BitmapFile = "dotmodel";

  if(getNewOne)  return new Dot_Model();
  return 0;
}
