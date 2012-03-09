/***************************************************************************
                          wirelabel.h  -  description
 ***************************************************************************/

#ifndef WIRELABEL_H
#define WIRELABEL_H

#include <qpainter.h>
#include <qstring.h>
#include <q3ptrlist.h>

#include "element_gui.h"
#include "viewpainter.h"

class Wire;
class Node;

class WireLabel : public Element_Qucs {
public:
  WireLabel(const QString& _Name=0, int _cx=0, int _cy=0,
            int _x1=0, int _y1=0, int _Type=isNodeLabel);
  ~WireLabel();

  void paintScheme(QPainter *p);
  void setCenter(int x, int y, bool relative=false);
  bool getSelected(int, int);
  void setName(const QString& Name_);

  Wire    *pWire;
  Node    *pNode;
  QString Name, initValue;

  void    paint(ViewPainter*);
  void    rotate();
  QString save();
  bool    load(const QString& s);
  bool    isHorizontal();
};

#endif
