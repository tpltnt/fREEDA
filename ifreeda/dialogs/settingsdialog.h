/***************************************************************************
                          settingsdialog.h  -  description
                             -------------------
    begin                : Mon Oct 20 2003
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

#ifndef SETTINGSDIALOG_H
#define SETTINGSDIALOG_H

#include <qdialog.h>
//Added by qt3to4:
#include <Q3VBoxLayout>

class QucsDoc;
class QLineEdit;
class QCheckBox;
class Q3VBoxLayout;
class QRegExpValidator;

/**
  *@author Michael Margraf
  */

class SettingsDialog : public QDialog  {
   Q_OBJECT
public:
  SettingsDialog(QucsDoc *d, QWidget *parent=0, const char *name=0);
  ~SettingsDialog();

private slots:
  void slotOK();
  void slotApply();

public:
  QucsDoc   *Doc;

  QLineEdit *Input_DataSet, *Input_DataDisplay;
  QLineEdit *Input_GridX, *Input_GridY;
  QCheckBox *Check_OpenDpl, *Check_GridOn;

  Q3VBoxLayout *all;
  QRegExpValidator *valExpr;
};

#endif
