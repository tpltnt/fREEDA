/***************************************************************************
                      guisettingsdialog.h  -  description
  ***************************************************************************/

#ifndef GUISETTINGSDIALOG_H
#define GUISETTINGSDIALOG_H

#include "gui.h"

#include <qdialog.h>
#include <qfont.h>
//Added by qt3to4:
#include <Q3VBoxLayout>

class QLineEdit;
class Q3VBoxLayout;
class QPushButton;
class QIntValidator;

/**
  *@author Michael Margraf
  */

class QucsSettingsDialog : public QDialog  {
   Q_OBJECT
public:
  QucsSettingsDialog(QucsApp *parent=0, const char *name=0);
  ~QucsSettingsDialog();

private slots:
  void slotOK();
  void slotApply();
  void slotFontDialog();
  void slotBGColorDialog();
  void slotDefaultValues();

public:
  QucsApp *App;

  QFont    Font;
  QPushButton *FontButton, *BGColorButton;
  QLineEdit   *undoNumEdit, *editorEdit;

  Q3VBoxLayout   *all;
  QIntValidator *val200;
};

#endif
