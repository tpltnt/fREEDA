/***************************************************************************
                          top.h  -  description
***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#ifndef DOT_TOP_H
#define DOT_TOP_H

#include "component.h"


class Dot_Top : public Component  {
public:
  Dot_Top();
  ~Dot_Top();
  Component* newOne();
  static Component* info(QString&, char* &, bool getNewOne=false);
};

#endif
