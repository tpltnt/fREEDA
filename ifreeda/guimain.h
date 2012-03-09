/***************************************************************************
                          guimain.h  -  description
 ***************************************************************************/

#ifndef GUIMAIN_H
#define GUIMAIN_H

#include <qfont.h>
#include <qcolor.h>
#include <q3ptrlist.h>

#include "wire.h"
#include "node.h"
#include "diagrams/diagram.h"

class QucsApp;
class Component;

// constants may be missing on windows systems
#include <math.h>
#ifndef M_PI
#define M_PI     3.1415926535897932384626433832795029
#endif

struct tQucsSettings {
  int x, y, dx, dy;    // position and size of main window
  QFont font;
  float largeFontSize;
  QColor BGColor;      // background color of view area

  unsigned int maxUndo;    // size of undo stack
  QString Editor;
  QString BinDir;
  QString BitmapDir;
  QString LangDir;
};

extern tQucsSettings QucsSettings;  // extern because nearly everywhere used
extern QFont savingFont;       // to remember which font to save in "qucsrc"

extern QucsApp *QucsMain;    // the Qucs application itself

bool saveApplSettings(QucsApp*);

QString complexRect(double, double, int Precision=3);
QString complexDeg (double, double, int Precision=3);
QString complexRad (double, double, int Precision=3);
QString StringNum  (double num, char form='g', int Precision=3);
void    str2num    (const QString&, double&, QString&, double&);

// just dummies for empty lists
extern Q3PtrList<Wire>      SymbolWires;
extern Q3PtrList<Node>      SymbolNodes;
extern Q3PtrList<Diagram>   SymbolDiags;
extern Q3PtrList<Component> SymbolComps;

#endif

