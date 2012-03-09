/***************************************************************************
                          element.cpp  -  description
                             -------------------
    begin                : Sat Sep 20 2003
    copyright            : (C) 2003 by Michael Margraf
    email                : michael.margraf@alumni.tu-berlin.de
 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#include "element_gui.h"

Element_Qucs::Element_Qucs()
{
  Type = isDummy;
}

Element_Qucs::~Element_Qucs()
{
}

void Element_Qucs::paintScheme(QPainter *)
{
}

void Element_Qucs::setCenter(int, int, bool)
{
}

void Element_Qucs::getCenter(int&, int&)
{
}
