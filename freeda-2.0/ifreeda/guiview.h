/***************************************************************************
                          guiview.h  -  description
 ***************************************************************************/

#ifndef GUIVIEW_H
#define GUIVIEW_H

#include <q3scrollview.h>
#include <q3ptrlist.h>
#include <qstring.h>
#include <qevent.h>
#include <qregexp.h>
#include <qvalidator.h>
#include <qlineedit.h>                                                                                                                                                                               
//Added by qt3to4:
#include <QWheelEvent>
#include <QMouseEvent>
#include <Q3PopupMenu>
#include <QScrollBar>
#include <QImage>

#include "components/component.h"
#include "guidoc.h"
#include "viewpainter.h"
#include "guimain.h"
#include <QRubberBand>

class Wire;
class QWidget;
class QPainter;
class Q3PopupMenu;


// **************************************************************************
// ******* This class provides an incomplete base for the application *******
// ******* view. (scematics, data displays etc.)                      *******
// **************************************************************************

class QucsView : public Q3ScrollView
{
  Q_OBJECT

public:
  QucsView(QWidget *parent=0);
  ~QucsView();

  float Zoom(float);
  bool   pasteElements();
  void   enlargeView(int, int, int, int);
  void   setPainter(QPainter*, QucsDoc*);
  void   eraseCross();
  void   editLabel(WireLabel*);

  Component *selComp;   // component selected in IconView
  Diagram   *selDiag;   // diagram selected in IconView
  Painting  *selPaint;  // painting selected in IconView

  bool    drawn;  // indicates whether the scheme element was drawn last time
  QString ProjName;
  QLineEdit *editText;  // for edit component properties

  Q3PtrList<Element_Qucs> movingElements;

  // menu appearing by right mouse button click on component
  Q3PopupMenu *ComponentMenu;


  // -------------------------------------------------------------------
  Q3PtrList<QucsDoc>  Docs; // document instances (schematics, data displays)

  //QImage main_image ;  // Added by Shivam. Painting is done in such a way that first viewport is drawn on scrollarea and then image is drawn on viewport.

//QPixmap main_image;

protected:
  void drawContents(QPainter*, int, int, int, int);
  void contentsMouseMoveEvent(QMouseEvent*);
  void contentsMousePressEvent(QMouseEvent*);
  void contentsMouseDoubleClickEvent(QMouseEvent*);
  void contentsMouseReleaseEvent(QMouseEvent*);
  void contentsWheelEvent(QWheelEvent*);

  bool ScrollUp(int);
  bool ScrollDown(int);
  bool ScrollLeft(int);
  bool ScrollRight(int);

protected slots:
  // changed by shivam
  // Added integer argument
  void slotScrollUp(int);
  void slotScrollDown(int);
  void slotScrollLeft(int);
  void slotScrollRight(int);
  void slotApplyCompText();

public slots:
  void slotCursorLeft();
  void slotCursorRight();
  void slotCursorUp();
  void slotCursorDown();

  void slotHideEdit();
  void slotEditElement();

public:
  void MouseDoNothing(QMouseEvent*);
  void MMoveSelect(QMouseEvent*);
  void MMoveComponent(QMouseEvent*);
  void MMoveDiagram(QMouseEvent*);
  void MMoveWire1(QMouseEvent*);
  void MMoveWire2(QMouseEvent*);
  void MMoveMoving(QMouseEvent*);
  void MMoveMoving2(QMouseEvent*);
  void MMovePaste(QMouseEvent*);
  void MMovePainting(QMouseEvent*);
  void MMoveDelete(QMouseEvent*);
  void MMoveLabel(QMouseEvent*);
  void MMoveMarker(QMouseEvent*);
  void MMoveMirrorY(QMouseEvent*);
  void MMoveMirrorX(QMouseEvent*);
  void MMoveRotate(QMouseEvent*);
  void MMoveActivate(QMouseEvent*);
  void MMoveOnGrid(QMouseEvent*);
  void MMoveResizePainting(QMouseEvent*);
  void MMoveMoveText(QMouseEvent*);
  void MMoveMoveTextB(QMouseEvent*);
  void MMoveZoomIn(QMouseEvent*);
  void (QucsView::*MouseMoveAction) (QMouseEvent*);// current mouse move method

  void MPressSelect(QMouseEvent*);
  void MPressDelete(QMouseEvent*);
  void MPressActivate(QMouseEvent*);
  void MPressMirrorX(QMouseEvent*);
  void MPressMirrorY(QMouseEvent*);
  void MPressRotate(QMouseEvent*);
  void MPressComponent(QMouseEvent*);
  void MPressDiagram(QMouseEvent*);
  void MPressLabel(QMouseEvent*);
  void MPressWire1(QMouseEvent*);
  void MPressWire2(QMouseEvent*);
  void MPressPainting(QMouseEvent*);
  void MPressMarker(QMouseEvent*);
  void MPressOnGrid(QMouseEvent*);
  void MPressMoveText(QMouseEvent*);
  void MPressZoomIn(QMouseEvent*);
  void (QucsView::*MousePressAction) (QMouseEvent*); // mouse press method

  void MDoubleClickSelect(QMouseEvent*);
  void MDoubleClickWire2(QMouseEvent*);
  void (QucsView::*MouseDoubleClickAction) (QMouseEvent*);

  void MReleaseSelect(QMouseEvent*);
  void MReleaseSelect2(QMouseEvent*);
  void MReleaseActivate(QMouseEvent*);
  void MReleaseMoving(QMouseEvent*);
  void MReleaseResizeDiagram(QMouseEvent*);
  void MReleasePaste(QMouseEvent*);
  void MReleaseResizePainting(QMouseEvent*);
  void MReleaseMoveText(QMouseEvent*);
  void MReleaseZoomIn(QMouseEvent*);
  void (QucsView::*MouseReleaseAction) (QMouseEvent*);

  void MovingElements();
  void endElementMoving();

private:
  void editElement(QMouseEvent*);
  void rightPressMenu(QMouseEvent*);
  void PressLabel(QMouseEvent*);

  int MAx1, MAy1,MAx2, MAy2, MAx3, MAy3;  // cache for mouse movements
  bool isMoveEqual;
  Element_Qucs *focusElement;
  QMouseEvent *focusMEvent;
  Wire *labeledWire;     // remember the wire whose label is moving
  ViewPainter Painter;

  QRegExp Expression;
  QRegExpValidator *Validator;

 // protected:
 // void paintEvent(QPaintEvent*);
};

#endif
