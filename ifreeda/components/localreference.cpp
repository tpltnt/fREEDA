/***************************************************************************
                          localreference.cpp  -  description
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

#include "localreference.h"


LocalReference::LocalReference()
{
  Description = QObject::tr("local (local reference potential)");

  Lines.append(new Line(  0,  0,  0, 10,QPen(Qt::darkBlue,2)));
  Lines.append(new Line(  0, 10,-10, 15,QPen(Qt::darkBlue,3)));
  Lines.append(new Line(-10, 15,  0, 20,QPen(Qt::darkBlue,3)));
  Lines.append(new Line(  0, 20, 10, 15,QPen(Qt::darkBlue,3)));
  Lines.append(new Line( 10, 15,  0, 10,QPen(Qt::darkBlue,3)));

  Ports.append(new Port(  0,  0));

  x1 = -12; y1 =  0;
  x2 =  12; y2 = 25;

  tx = 0;
  ty = 0;
  Model = "LRT";
  Name  = "";
}

LocalReference::~LocalReference()
{
}

Component* LocalReference::newOne()
{
  return new LocalReference();
}

Component* LocalReference::info(QString& Name, char* &BitmapFile, bool getNewOne)
{
  Name = QObject::tr("Local Reference");
  BitmapFile = "localreference";

  if(getNewOne)  return new LocalReference();
  return 0;
}

