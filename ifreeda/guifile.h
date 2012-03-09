/***************************************************************************
                          guifile.h  -  description
 ***************************************************************************/

#ifndef GUIFILE_H
#define GUIFILE_H


//Added by qt3to4:
#include <Q3TextStream>
#include <Q3PtrList>

#include "wire.h"
#include "diagrams/diagram.h"
#include "paintings/painting.h"
#include "components/component.h"


class Q3TextEdit;
class Q3Process;
class QucsDoc;

class QucsFile {
public:
  QucsFile(QucsDoc*);
  ~QucsFile();

  int   save();
  bool  load();

  QString createClipboardFile();
  bool    pasteFromClipboard(Q3TextStream*, Q3PtrList<Element_Qucs>*);
  QString createUndoString(char);
  bool    rebuild(QString *);
  QString createSymbolUndoString(char);
  bool    rebuildSymbol(QString *);

  bool  createSubNetlist(Q3TextStream*, int&, QStringList&, Q3TextEdit*);
  bool  prepareNetlist(Q3TextStream&, QStringList&, Q3TextEdit*);
  void  createNetlist(Q3TextStream&);


private:
  bool  loadProperties(Q3TextStream*);
  void  simpleInsertComponent(Component*);
  bool  loadComponents(Q3TextStream*, Q3PtrList<Component> *List=0);
  void  simpleInsertWire(Wire*);
  bool  loadWires(Q3TextStream*, Q3PtrList<Element_Qucs> *List=0);
  bool  loadDiagrams(Q3TextStream*, Q3PtrList<Diagram>*);
  bool  loadPaintings(Q3TextStream*, Q3PtrList<Painting>*);
  bool  loadIntoNothing(Q3TextStream*);

  bool  giveNodeNames(Q3TextStream*, int&, QStringList&, Q3TextEdit*);

  QucsDoc  *Doc;
  Q3PtrList<Wire>      *Wires;
  Q3PtrList<Node>      *Nodes;
  Q3PtrList<Component> *Comps;
  Q3PtrList<Diagram>   *Diags;
  Q3PtrList<Painting>  *Paints, *SymbolPaints;

};

#endif
