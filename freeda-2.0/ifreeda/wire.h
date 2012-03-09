/***************************************************************************
                          wire.h  -  description
 ***************************************************************************/

#ifndef WIRE_H
#define WIRE_H

#include <qstring.h>
#include <q3ptrlist.h>
#include <qpainter.h>


#include "viewpainter.h"
#include "element_gui.h"
#include "components/component.h"    // because of struct Port
#include "wirelabel.h"


class Wire : public Element_Qucs {
public:
  Wire(int _x1=0, int _y1=0, int _x2=0, int _y2=0, Node *n1=0, Node *n2=0);
  ~Wire();

  void paint(ViewPainter*);
  void paintScheme(QPainter*);
  void setCenter(int, int, bool relative=false);
  void getCenter(int&, int&);
  bool getSelected(int, int);
  void setName(const QString&, const QString&, int delta_=0, int x_=0, int y_=0);

  Node      *Port1, *Port2;
  WireLabel *Label;

  void    rotate();
  QString save();
  bool    load(const QString&);
  bool    isHorizontal();
};

#endif
