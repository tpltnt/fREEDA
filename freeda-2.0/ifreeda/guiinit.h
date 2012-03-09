/***************************************************************************
                          guiinit.h  -  description
 ***************************************************************************/

#ifndef GUIINIT_H
#define GUIINIT_H

#include <qaction.h>
#include <q3accel.h>
#include <qobject.h>
#include <q3toolbar.h>
#include <q3popupmenu.h>

class QucsApp;
class QucsActions;


/**
  *@author Michael Margraf
  */

class QucsInit : public QObject {
  Q_OBJECT
public:
  QucsInit();
  ~QucsInit();

  void perform(QucsApp *p_);

public slots:
  void slotViewToolBar(bool toggle);    // toggle the toolbar
  void slotViewStatusBar(bool toggle);  // toggle the statusbar

  void slotHelpAbout();       // shows an about dialog
  void slotHelpAboutQt();     // shows the standard about dialog for Qt

private:
  void initActions();    // initializes all QActions of the application
  void initMenuBar();    // creates the menu_bar and inserts the menuitems
  void initToolBar();    // creates the toolbars
  void initStatusBar();  // setup the statusbar

  QucsApp *App;   // the application that called this instance
  QucsActions *Acts;

  QAction *helpAboutApp, *helpAboutQt, *viewToolBar, *viewStatusBar;

  // menus contain the items of their menubar
  Q3PopupMenu *fileMenu, *editMenu, *insMenu, *projMenu, *simMenu, *viewMenu;
  Q3PopupMenu *helpMenu, *alignMenu, *toolMenu;

  Q3ToolBar *fileToolbar, *editToolbar, *viewToolbar, *workToolbar;

  Q3Accel *mainAccel;     // to set more than one key to one action
};

#endif
