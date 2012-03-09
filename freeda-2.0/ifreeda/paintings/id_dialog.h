/***************************************************************************
                          id_dialog.h  -  description
                             -------------------
    begin                : Sat Oct 16 2004
    copyright            : (C) 2004 by Michael Margraf
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

#ifndef ID_DIALOG_H
#define ID_DIALOG_H

#include <qdialog.h>
#include <qregexp.h>
//Added by qt3to4:
#include <Q3VBoxLayout>


class QLineEdit;
class Q3VBoxLayout;
class QRegExp;
class QRegExpValidator;

/**
  *@author Michael Margraf
  */

class ID_Dialog : public QDialog  {
Q_OBJECT
public:
  ID_Dialog(QWidget *parent=0);
  ~ID_Dialog();

public:
  QLineEdit   *Prefix;
  Q3VBoxLayout *v;

  QRegExp *rx;
  QRegExpValidator *Validator;
};

#endif
