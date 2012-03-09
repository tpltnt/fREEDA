/***************************************************************************
                          guiactions.h  -  description
 ***************************************************************************/

#ifndef GUIACTIONS_H
#define GUIACTIONS_H

#include <qaction.h>
//Added by qt3to4:
#include <QMouseEvent>


class QucsApp;
class QucsDoc;
class QucsView;



typedef bool (QucsDoc::*pToggleFunc) (); // pointer to toggle action
typedef void (QucsView::*pMouseFunc) (QMouseEvent*);


class QucsActions : public QObject {
  Q_OBJECT
public:
  QucsActions();
  ~QucsActions();

  void init(QucsApp *p_);
  void editFile(const QString&);

  QAction *insWire, *insLabel, *insGround, *insPort, *insPlot, *insDot_Model, *insDot_Top,
          *insDot_Bottom, *insEquation, *magPlus, *editRotate, *editMirror, *editMirrorY,
	  *editPaste, *select, *editActivate, *wire, *editDelete, *setMarker, *onGrid,
	  *moveText, *helpIndex, *helpGetStart, *callEditor, *callFilter, *callLine,
	  *showMsg, *showNet, *alignTop, *alignBottom, *alignLeft, *alignRight, *distrHor,
	  *distrVert, *selectAll;

public slots:
  void slotEditRotate(bool);  // rotate the selected items
  void slotEditMirrorX(bool); // mirror the selected items about X axis
  void slotEditMirrorY(bool); // mirror the selected items about Y axis
  void slotEditPaste(bool);   // paste the clipboard into the document
  void slotEditDelete(bool);  // delete the selected items
  void slotInsertEquation(bool);
  void slotInsertGround(bool);
  void slotInsertPort(bool);
  void slotInsertPlot(bool);
  void slotInsertDot_Model(bool);
  void slotInsertDot_Top(bool);
  void slotInsertDot_Bottom(bool);
  void slotSetWire(bool);
  void slotSelect(bool);
  void slotEditActivate(bool);
  void slotInsertLabel(bool);
  void slotSetMarker(bool);
  void slotOnGrid(bool);      // set selected elements on grid
  void slotMoveText(bool);    // move property text of components
  void slotZoomIn(bool);
  void slotEditUndo();    // makes the last operation undone
  void slotEditRedo();    // makes the last undo undone
  void slotAlignTop();    // align selected elements with respect to top
  void slotAlignBottom(); // align selected elements with respect to bottom
  void slotAlignLeft();   // align selected elements with respect to left
  void slotAlignRight();  // align selected elements with respect to right
  void slotDistribHoriz();// distribute horizontally selected elements
  void slotDistribVert(); // distribute vertically selected elements
  void slotSelectAll();
  void slotShowLastMsg();
  void slotShowLastNetlist();
  void slotCallEditor();
  void slotCallFilter();
  void slotCallLine();
  void slotHelpIndex();       // shows a HTML docu: Help Index
  void slotGettingStarted();  // shows a HTML docu: Getting started

private:
  void showHTML(const QString&);
  bool performToggleAction(bool, QAction*, pToggleFunc, pMouseFunc, pMouseFunc);

  // copies of variables in QucsApps
  QucsApp  *App;
  QucsView *view;
};

#endif
