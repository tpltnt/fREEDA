/***************************************************************************
                          guidoc.h  -  description
 ***************************************************************************/

#ifndef GUIDOC_H
#define GUIDOC_H



#include <qstring.h>
#include <qpainter.h>
#include <q3ptrlist.h>
#include <qfile.h>
#include <qtabbar.h>
//Added by qt3to4:
#include <Q3TextStream>

#include "wire.h"
#include "diagrams/diagram.h"
#include "paintings/painting.h"
#include "components/component.h"
#include "guifile.h"


class QucsApp;


class QucsDoc {
public: 
  QucsDoc(QucsApp*, const QString&);
  ~QucsDoc();

  void setName(const QString&);
  void setChanged(bool, bool fillStack=false, char Op='*');

  void paint(ViewPainter*);
  void paintGrid(ViewPainter*, int, int, int, int);
  void print(QPainter*, bool);

  int   insertWireNode1(Wire*);
  bool  connectHWires1(Wire*);
  bool  connectVWires1(Wire*);
  int   insertWireNode2(Wire*);
  bool  connectHWires2(Wire*);
  bool  connectVWires2(Wire*);
  int   insertWire(Wire*);

  Node* insertNode(int, int, Element_Qucs*);
  void  insertComponentNodes(Component*);
  void  insertRawComponent(Component*, bool num=false);
  void  insertComponent(Component*);
  void  insertNodeLabel(WireLabel*);

  Component* searchSelSubcircuit();
  void       sizeOfAll(int&, int&, int&, int&);
  Component* selectedComponent(int, int);
  Node*      selectedNode(int, int);
  Wire*      selectedWire(int, int);
  Painting*  selectedPainting(int, int);
  void       selectWireLine(Element_Qucs*, Node*, bool);
  Element_Qucs*   selectElement(int, int, bool, int *index=0);
  int        selectElements(int, int, int, int, bool);
  void       deselectElements(Element_Qucs*);
  bool  activateComponent(int, int);
  bool  activateComponents();
  void  activateComps(int, int, int, int);
  void  NewMovingWires(Q3PtrList<Element_Qucs>*, Node*);
  void  copySelectedElements(Q3PtrList<Element_Qucs>*);

  void  setComponentNumber(Component*);
  void  oneLabel(Node*);
  int   placeNodeLabel(WireLabel*);
  Element_Qucs* getWireLabel(Node*);
  void  setCompPorts(Component*);
  void  copyComponents(int&, int&, int&, int&);
  void  copyComponents2(int&, int&, int&, int&);
  void  copyWires(int&, int&, int&, int&);
  void  copyLabels(int&, int&, int&, int&);
  void  copyPaintings(int&, int&, int&, int&);
  bool  copyComps2WiresPaints(int&, int&, int&, int&);
  bool  rotateElements();
  bool  mirrorXComponents();
  bool  mirrorYComponents();
  bool  oneTwoWires(Node*);
  void  deleteComp(Component*);
  void  deleteWire(Wire*);
  bool  deleteElements();
  Marker*  setMarker(int, int);
  bool  MarkerLeftRight(bool);
  bool  MarkerUpDown(bool);
  int   copyElements(int&, int&, int&, int&);
  bool  aligning(int);
  bool  distribHoriz();
  bool  distribVert();
  bool  elementsOnGrid();
  void  switchPaintMode();

  QString copySelected(bool);
  bool    paste(Q3TextStream*, Q3PtrList<Element_Qucs>*);
  bool    load();
  int     save();
  int     adjustPortNumbers();
  bool    undo();
  bool    redo();

  void    reloadGraphs();
  void    setOnGrid(int&, int&);
  Component* selectCompText(int, int, int&, int&);


  QucsFile  File;   // class to perform  load, save, copy, paste

  QString DocName;
  bool    DocChanged;

  QucsApp *App;
  //QTab    *Tab;
  QTabBar *Bar;

  Q3PtrList<Wire>      *Wires,  DocWires;
  Q3PtrList<Node>      *Nodes,  DocNodes;
  Q3PtrList<Diagram>   *Diags,  DocDiags;
  Q3PtrList<Painting>  *Paints, DocPaints;
  Q3PtrList<Component> *Comps,  DocComps;

  Q3PtrList<Painting>  SymbolPaints;  // symbol definition for subcircuit
  bool  symbolMode;  // true if in symbol painting mode

  bool    SimOpenDpl;  // open data display after simulation ?
  bool    GenNetlist;  // Only generate the netlist do not simulate
  QString DataSet;     // name of the default dataset
  QString DataDisplay; // name of the default data display
  int     GridX, GridY;
  bool    GridOn;

  float  Scale;
  int PosX, PosY; // upper left corner of visible area (only for remembering during seeing another document)
  int ViewX1, ViewY1, ViewX2, ViewY2;  // size of the document area
  int UsedX1, UsedY1, UsedX2, UsedY2;  // document area used by elements

  // Two of those data sets are needed for Schematic and for symbol.
  // Which one is in "tmp..." depends on "symbolMode".
  float  tmpScale;
  int tmpPosX, tmpPosY;
  int tmpViewX1, tmpViewY1, tmpViewX2, tmpViewY2;
  int tmpUsedX1, tmpUsedY1, tmpUsedX2, tmpUsedY2;
  int tab_new;  
  Q3PtrList<Element_Qucs> ElementCache;
  Q3PtrList<QString> UndoStack;
  Q3PtrList<QString> UndoSymbol;    // undo stack for circuit symbol
};


#endif
