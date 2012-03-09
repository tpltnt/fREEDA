/***************************************************************************
                          qucs.cpp  -  description
                             -------------------
    begin                : Thu Aug 28 18:17:41 CEST 2003
    copyright            : (C) 2003, 2004 by Michael Margraf
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
#include<qlocale.h>
#include<qabstractitemmodel.h>
// Above lines added by Shivam to correct header file errors

#include <q3accel.h>
#include <qimage.h>
#include <qsplitter.h>
#include <q3vbox.h>
#include <q3hbox.h>
#include <qlabel.h>
#include <qmessagebox.h>
#include <qdir.h>
#include <qpainter.h>
#include <q3filedialog.h>
#include <qinputdialog.h>
#include <qclipboard.h>
#include <qfont.h>
#include <q3textedit.h>
#include <qcheckbox.h>
//Added by qt3to4:
#include <QPixmap>
#include <QCloseEvent>
#include <Q3TextStream>
#include <Q3PopupMenu>
#include <string.h>
using namespace std;
#include <config.h>
#include <limits.h>

#include "guimain.h"
#include "gui.h"
#include "guiview.h"
#include "wire.h"
#include "components/components.h"
#include "paintings/paintings.h"
#include "diagrams/diagrams.h"
#include "dialogs/messagebox.h"
#include "dialogs/newprojdialog.h"
#include "dialogs/settingsdialog.h"
#include "dialogs/guisettingsdialog.h"
#include "dialogs/simmessage.h"

#include <iostream>


#define  COMBO_fREEDA_lumped    0   
#define  COMBO_fREEDA_nonlinear 1
#define  COMBO_fREEDA_sources   2
#define  COMBO_fREEDA_logic     3
#define  COMBO_fREEDA_thermal   4
#define  COMBO_fREEDA_other     5
#define  COMBO_fanalysis        6
#define  COMBO_Diagrams         7
#define  COMBO_File             8
//***Old qucs categories dont display *********
#define  COMBO_passive   9
#define  COMBO_Sources   10
#define  COMBO_TLines    11
#define  COMBO_nonlinear 12
#define  COMBO_Sims      13
//*********************************************
#define  COMBO_Paints    14  // must be the last one

QucsApp::QucsApp()
{
  setCaption("ifREEDA " PACKAGE_VERSION);

  QucsFileFilter = tr("Schematic")+" (*.sch);;"+
                   tr("Data Display")+" (*.dpl);;"+
		   tr("Qucs Documents")+" (*.sch *.dpl);;"+
		   tr("Any File")+" (*)";
  QucsWorkDir.setPath(QDir::currentDirPath()+QDir::convertSeparators("/projects"));
  QucsHomeDir.setPath(QDir::currentDirPath()+QDir::convertSeparators("/projects"));

  move  (QucsSettings.x,  QucsSettings.y);
  resize(QucsSettings.dx, QucsSettings.dy);
//  resize(maximumSize());

  initView();
  Init.perform(this);  // initializes toolbar, statusbar, menubar, actions
  Acts.init(this);
  initCursorMenu();


  // default settings of the printer
  Printer = new QPrinter(QPrinter::PrinterResolution);
  Printer->setOrientation(QPrinter::Landscape);
  Printer->setColorMode(QPrinter::Color);

  Acts.select->setOn(true);  // switch on the 'select' action
  HierarchyHistory.setAutoDelete(true);

  // creates a document called "untitled"
  view->Docs.append(new QucsDoc(this, ""));
  //cout << "Printing Untitled pointer" << view->Docs.current() << endl;

  // load documents given as command line arguments
  for(int z=1; z<qApp->argc(); z++)
    if(*(qApp->argv()[z]) != '-')
      gotoPage(qApp->argv()[z]);
}

QucsApp::~QucsApp()
{
  delete Printer;
}


// #########################################################################
// ##########                                                     ##########
// ##########       Creates the view area (QTabWidget etc.)       ##########
// ##########                                                     ##########
// #########################################################################
void QucsApp::initView()
{
  // set application icon
  setIcon (QPixmap(QucsSettings.BitmapDir + "big.qucs.xpm"));
  Q3VBox *all = new Q3VBox(this);   // only to fill the entire view area
  QSplitter *Hsplit = new QSplitter(Qt::Horizontal, all);

  TabView  = new QTabWidget(Hsplit);    // tabs on the left side
  Q3VBox *WorkGroup = new Q3VBox(Hsplit);
  WorkView = new QTabBar(WorkGroup);    // tab on the right side
  view = new QucsView(WorkGroup);    // work area with documents

  connect(WorkView, SIGNAL(selected(int)), SLOT(slotChangeView(int)));

  Hsplit->setResizeMode(TabView, QSplitter::KeepSize);
  setCentralWidget(all);

  // ----------------------------------------------------------
  // "Project Tab" of the left QTabWidget
  Q3VBox *ProjGroup = new Q3VBox(this);
  Q3HBox *ProjButts = new Q3HBox(ProjGroup);
  QPushButton *ProjNew   = new QPushButton(tr("New"),ProjButts);
  connect(ProjNew, SIGNAL(clicked()), SLOT(slotProjNewButt()));
  QPushButton *ProjOpen  = new QPushButton(tr("Open"),ProjButts);
  connect(ProjOpen, SIGNAL(clicked()), SLOT(slotProjOpenButt()));
  QPushButton *ProjDel   = new QPushButton(tr("Delete"),ProjButts);
  connect(ProjDel, SIGNAL(clicked()), SLOT(slotProjDelButt()));

  Projects = new Q3ListBox(ProjGroup);
  TabView->addTab(ProjGroup, tr("Projects"));
  TabView->setTabToolTip(TabView->page(0),
			 tr("content of the project directory"));

  connect(Projects, SIGNAL(doubleClicked(Q3ListBoxItem*)),
		    SLOT(slotOpenProject(Q3ListBoxItem*)));

  // ----------------------------------------------------------
  // "Content Tab" of the left QTabWidget
  Content = new Q3ListView(this);
  Content->setRootIsDecorated(true); // open/close decoration for root items
  Content->setSorting(-1);    // no sorting
  Content->addColumn(tr("Content of"));
  Content->addColumn(tr("Note"));
  Content->setColumnWidthMode(0,Q3ListView::Manual);
  Content->setColumnWidth(0, 150);

  ConDatasets   = new Q3ListViewItem(Content, tr("Datasets"));
  ConDisplays   = new Q3ListViewItem(Content, tr("Data Displays"));
  ConSchematics = new Q3ListViewItem(Content, tr("Schematics"));
  TabView->addTab(Content,tr("Content"));
  TabView->setTabToolTip(TabView->page(1), tr("content of the open project"));

// QT 3.2
//  connect(Content, SIGNAL(doubleClicked(QListViewItem*, const QPoint &,int)), SLOT(slotOpenContent(QListViewItem*, const QPoint &,int)));
  connect(Content, SIGNAL(doubleClicked(Q3ListViewItem*)),
		   SLOT(slotOpenContent(Q3ListViewItem*)));
  connect(Content, SIGNAL(clicked(Q3ListViewItem*)),
		   SLOT(slotSelectSubcircuit(Q3ListViewItem*)));

  // ----------------------------------------------------------
  // "Component Tab" of the left QTabWidget
  Q3VBox *CompGroup  = new Q3VBox(this);
  CompChoose = new QComboBox(CompGroup);
  CompComps  = new Q3IconView(CompGroup);
  TabView->addTab(CompGroup,tr("Components"));
  TabView->setTabToolTip(TabView->page(2), tr("components and diagrams"));
  fillComboBox(true);

  slotSetCompView(0);
  connect(CompChoose, SIGNAL(activated(int)), SLOT(slotSetCompView(int)));
  connect(CompComps, SIGNAL(clicked(Q3IconViewItem*)),
		     SLOT(slotSelectComponent(Q3IconViewItem*)));


  // ---------------------------------------------------------------------
  readProjects();   // reads all projects and inserts them into the ListBox
}

// ----------------------------------------------------------
void QucsApp::fillComboBox(bool setAll)
{
  CompChoose->clear();
  if(setAll) {
    CompChoose->insertItem(tr("fREEDA Lumped"));
    CompChoose->insertItem(tr("fREEDA Nonlinear"));
    CompChoose->insertItem(tr("fREEDA Sources"));
    CompChoose->insertItem(tr("fREEDA Logic"));
    CompChoose->insertItem(tr("fREEDA Thermal"));
    CompChoose->insertItem(tr("fREEDA Other"));
    CompChoose->insertItem(tr("fREEDA Analysis"));

    //CompChoose->insertItem(tr("lumped components"));   //QUCS-EDIT
    //CompChoose->insertItem(tr("sources"));             //QUCS-EDIT
   // CompChoose->insertItem(tr("transmission lines"));  //QUCS-EDIT
   // CompChoose->insertItem(tr("nonlinear components"));//QUCS-EDIT
   // CompChoose->insertItem(tr("simulations"));         //QUCS-EDIT
    CompChoose->insertItem(tr("diagrams"));
    CompChoose->insertItem(tr("file data"));

  }
  CompChoose->insertItem(tr("paintings"));
 
}

// ----------------------------------------------------------
// Menu that appears if right mouse button is pressed on a file in the
// "Content" ListView.
void QucsApp::initCursorMenu()
{
  ContentMenu = new Q3PopupMenu(Content);
  ContentMenu->insertItem(tr("Open"), this, SLOT(slotCMenuOpen()));
  ContentMenu->insertItem(tr("Rename"), this, SLOT(slotCMenuRename()));
  ContentMenu->insertItem(tr("Delete"), this, SLOT(slotCMenuDelete()));
  ContentMenu->insertItem(tr("Delete Group"), this, SLOT(slotCMenuDelGroup()));

  connect(Content,
	  SIGNAL(contextMenuRequested(Q3ListViewItem*, const QPoint&, int)),
	  SLOT(slotShowContentMenu(Q3ListViewItem*, const QPoint&, int)));
}

// ----------------------------------------------------------
// Is called if right mouse button is pressed on a file in the "Content" ListView.
// It shows a menu.
void QucsApp::slotShowContentMenu(Q3ListViewItem *item, const QPoint& point, int)
{
  if(item)
    if(item->parent() != 0)   // no component, but item "schematic", ...
      ContentMenu->popup(point);
}

// ----------------------------------------------------------
// for menu that appears if right mouse button is pressed on a file in the "Content" ListView.
void QucsApp::slotCMenuOpen()
{
  Q3ListViewItem *i = Content->selectedItem();
  if(i == 0) return;

  slotOpenContent(i);
}

// ----------------------------------------------------------
// for menu that appears if right mouse button is pressed on a file in the "Content" ListView.
void QucsApp::slotCMenuRename()
{
  Q3ListViewItem *i = Content->selectedItem();
  if(i == 0) return;

  QucsDoc *dc = view->Docs.current();
  // search, if file is open
  for(QucsDoc *d = view->Docs.first(); d!=0; d = view->Docs.next()) {
    if(d->DocName == QucsWorkDir.filePath(i->text(0))) {
      QMessageBox::critical(this, tr("Error"),
				  tr("Cannot rename an open file!"));
      view->Docs.findRef(dc);
      return;
    }
  }
  view->Docs.findRef(dc);

  QString Name   = i->text(0);
  QString Suffix = Name.right(4);    // remember suffix

  bool ok;
  QString s = QInputDialog::getText(tr("Rename file"), tr("Enter new name:"),
		QLineEdit::Normal, Name.left(Name.length()-4), &ok, this);
  if(!ok) return;
  if(s.isEmpty()) return;

  QDir file(QucsWorkDir.path());
  if(!file.rename(Name, s+Suffix)) {
    QMessageBox::critical(this, tr("Error"), tr("Cannot rename file: ")+Name);
    return;
  }
  i->setText(0, s+Suffix);
}

// ----------------------------------------------------------
// for menu that appears if right mouse button is pressed on a file in
// the "Content" ListView.
void QucsApp::slotCMenuDelete()
{
  Q3ListViewItem *item = Content->selectedItem();
  if(item == 0) return;
  QString FileName = QucsWorkDir.filePath(item->text(0));

  QucsDoc *dc = view->Docs.current();
  // search, if file is open
  for(QucsDoc *d = view->Docs.first(); d!=0; d = view->Docs.next()) {
    if(d->DocName == FileName) {
      QMessageBox::critical(this, tr("Error"),
                   tr("Cannot delete an open file!"));
      view->Docs.findRef(dc);
      return;
    }
  }
  view->Docs.findRef(dc);  // back to the real current document

  int n = QMessageBox::warning(this, tr("Warning"),
                 tr("This will delete the file permanently! Continue ?"),
                 tr("No"), tr("Yes"));
  if(n != 1) return;

  if(!QFile::remove(FileName)) {
    QMessageBox::critical(this, tr("Error"),
		tr("Cannot delete schematic: ")+item->text(0));
    return;
  }
//  if(ConSchematics == item->parent())
//    ConSchematics->takeItem(item);
//  else if(ConDisplays == item->parent())
//    ConDisplays->takeItem(item);
//  else
//    ConDatasets->takeItem(item);

  item->parent()->takeItem(item);
  delete item;
}

// ----------------------------------------------------------
// for menu that appears if right mouse button is pressed on a file in
// the "Content" ListView
// Deletes all files with that name (with suffix "sch", "dpl", "dat").
void QucsApp::slotCMenuDelGroup()
{
  Q3ListViewItem *item = Content->selectedItem();
  if(item == 0) return;
  QString s = item->text(0);
  s = s.left(s.length()-4);  // cut off suffix from file name

  QString NameSCH = QucsWorkDir.filePath(s+".sch");
  QString NameDPL = QucsWorkDir.filePath(s+".dpl");
  QString NameDAT = QucsWorkDir.filePath(s+".dat");

  QucsDoc *dc = view->Docs.current();
  // search, if files are open
  for(QucsDoc *d = view->Docs.first(); d!=0; d = view->Docs.next()) {
    if(d->DocName == NameSCH) {
      QMessageBox::critical(this, tr("Error"),
                   tr("Cannot delete the open file: ")+NameSCH);
      view->Docs.findRef(dc);
      return;
    }
    if(d->DocName == NameDPL) {
      QMessageBox::critical(this, tr("Error"),
                   tr("Cannot delete the open file: ")+NameDPL);
      view->Docs.findRef(dc);
      return;
    }
  }
  view->Docs.findRef(dc);  // back to the real current document

  bool SCH_exists = QFile::exists(NameSCH);
  bool DPL_exists = QFile::exists(NameDPL);
  bool DAT_exists = QFile::exists(NameDAT);

  QString Str;
  if(SCH_exists)  Str += s+".sch\n";
  if(DPL_exists)  Str += s+".dpl\n";
  if(DAT_exists)  Str += s+".dat\n";

  int n = QMessageBox::warning(this, tr("Warning"),
	  tr("This will delete the files\n")+Str+tr("permanently! Continue ?"),
	  tr("No"), tr("Yes"));
  if(n != 1) return;

  // remove files .................
  if(SCH_exists)
    if(!QFile::remove(NameSCH)) {
      QMessageBox::critical(this, tr("Error"),
		tr("Cannot delete schematic: ")+s+".sch");
      return;
    }
  if(DPL_exists)
    if(!QFile::remove(NameDPL)) {
      QMessageBox::critical(this, tr("Error"),
		tr("Cannot delete schematic: ")+s+".dpl");
      return;
    }
  if(DAT_exists)
    if(!QFile::remove(NameDAT)) {
      QMessageBox::critical(this, tr("Error"),
		tr("Cannot delete schematic: ")+s+".dat");
      return;
    }

  // remove items from listview ........
  if(SCH_exists) {
    item = Content->findItem(s+".sch", 0);
    if(item) {
      item->parent()->takeItem(item);
      delete item;
    }
  }
  if(DPL_exists) {
    item = Content->findItem(s+".dpl", 0);
    if(item) {
      item->parent()->takeItem(item);
      delete item;
    }
  }
  if(DAT_exists) {
    item = Content->findItem(s+".dat", 0);
    if(item) {
      item->parent()->takeItem(item);
      delete item;
    }
  }
}

// ----------------------------------------------------------
// Checks situation and reads all existing Qucs projects
void QucsApp::readProjects()
{
  QDir ProjDir(QucsHomeDir.path());

  if(!ProjDir.cd(QucsHomeDir.path())) // work directory exists ?
    if(!ProjDir.mkdir(QucsHomeDir.path())) { // no, then create it
      QMessageBox::warning(this, tr("Warning"),
			tr("Cannot create work directory!"));
      return;
    }
  ProjDir.cd(QucsHomeDir.path());

  // get all directories
  QStringList PrDirs = ProjDir.entryList("*", QDir::Dirs, QDir::Name);
  PrDirs.pop_front(); // delete "." from list
  PrDirs.pop_front(); // delete ".." from list

  Projects->clear();
  QStringList::iterator it;
  // inserts all project directories
  for(it = PrDirs.begin(); it != PrDirs.end(); it++)
     if ((*it).right(4) == "_prj") {   // project directories end with "_prj"
      (*it) = (*it).left((*it).length()-4);// remove "_prj" from project name
      Projects->insertItem(*it);
     }
}

// #########################################################################
bool QucsApp::closeAllFiles()
{
  int  Result = 0;
  bool notForAll = true;
  CloseMessageBox *m = new CloseMessageBox(tr("Closing ifREEDA document"),
	tr("This document contains unsaved changes!\n"
	   "Do you want to save the changes before closing?"),this);

  // close all files and ask to save changed ones
  for(QucsDoc *ptr = view->Docs.first(); ptr != 0; ) {
    nextDocument(false);
    if(ptr->DocChanged) {
      if(notForAll)  Result = m->exec();
      switch(Result) {
	case 1: if(!saveCurrentFile())  return false;  // save button
		break;
	case 2: Result = 1;         // save all button
		notForAll = false;
		if(!saveCurrentFile())  return false;
		break;
	case 4: Result = 3;         // discard all button
		notForAll = false;
		break;
	case 5: return false;       // cancel button
      }
    }

    view->Docs.remove();
    ptr = view->Docs.current();
  }

  delete m;
  switchEditMode(true);   // set schematic edit mode
  return true;
}

// ########################################################################
// Switches to the next document. Is called when closing a document.
void QucsApp::nextDocument(bool loadDiagrams)
{
  // Added by Shivam
     //WorkView->setCurrentTab(WorkView->tabAt(QPoint(view->Docs.at(),0)));
       WorkView->setCurrentIndex(view->Docs.at()); // setCurrentIndex(index) sets currentIndex property to index which reflects tab bar's visible tab   	

  // make new document the current
  //WorkView->setCurrentTab(WorkView->tabAt(view->Docs.at()));

  QucsDoc *d = view->Docs.current();
  view->resizeContents (int(d->Scale*double(d->ViewX2-d->ViewX1)),
			int(d->Scale*double(d->ViewY2-d->ViewY1)));
  view->setContentsPos(d->PosX, d->PosY);   // set view area

  if(loadDiagrams)
    view->Docs.current()->reloadGraphs();  // load recent simulation data
  view->viewport()->repaint();
  view->drawn = false;

  QString *ps = d->UndoStack.current();
  if(ps != d->UndoStack.getFirst())  undo->setEnabled(true);
  else  undo->setEnabled(false);
  if(ps != d->UndoStack.getLast())  redo->setEnabled(true);
  else  redo->setEnabled(false);

  HierarchyHistory.clear();   // no subcircuit history
  popH->setEnabled(false);
}

// ########################################################################
// Is called when another document is selected via the TabBar.
void QucsApp::slotChangeView(int id)
{
    view->editText->setHidden(true); // disable text edit of component property
    //cout << "Shivam" << endl; 
    QucsDoc *d = view->Docs.current();
    if(d !=0)
    {
    d->PosX = view->contentsX();    // save position for old document
    d->PosY = view->contentsY();

  // Changed by Shivam : In Qt4.3+  index, id and position index means the same. In Qt 3.3 id is returned after addTab command. 
     d = view->Docs.at(id);
 //  d = view->Docs.at(WorkView->indexOf(id));   // new current document

  view->resizeContents(int(d->Scale*double(d->ViewX2-d->ViewX1)),
                       int(d->Scale*double(d->ViewY2-d->ViewY1)));
  view->setContentsPos(d->PosX, d->PosY);

  // which mode: schematic or symbol editor ?
  if((CompChoose->count() > 1) == d->symbolMode)
    changeSchematicSymbolMode(d);

  view->Docs.current()->reloadGraphs();  // load recent simulation data
  view->viewport()->repaint();
  view->drawn = false;


  QString *ps = d->UndoStack.current();
  if(ps != d->UndoStack.getFirst())  undo->setEnabled(true);
  else  undo->setEnabled(false);
  if(ps != d->UndoStack.getLast())  redo->setEnabled(true);
  else  redo->setEnabled(false);

  if(!HierarchyHistory.isEmpty())
    if(*(HierarchyHistory.getLast()) != "*") {
      HierarchyHistory.clear();   // no subcircuit history anymore
      popH->setEnabled(false);
    }
    }
}

// #######################################################################
// Changes to the document "Name" if possible.
bool QucsApp::gotoPage(const QString& Name)
{
  QucsDoc *d;
  int cNo, No = view->Docs.at();

  // search, if page is already loaded
  for(d = view->Docs.first(); d!=0; d = view->Docs.next())
    if(d->DocName == Name) break;

  if(d != 0) {   // open page found ?
    view->Docs.current()->reloadGraphs();   // load recent simulation data
    // make new document the current
    //WorkView->setCurrentTab(WorkView->tabAt(view->Docs.at()));
    // Added by Shivam
     WorkView->setCurrentIndex(view->Docs.at());
    return true;
  }

  d = new QucsDoc(this, Name);
  view->Docs.append(d);   // create new page

  if(!d->load()) {    // load document if possible
    view->Docs.remove();
    view->Docs.at(No);
    view->viewport()->repaint();
    view->drawn = false;
    return false;
  }

  cNo = view->Docs.at();
  view->Docs.at(No);  // back to the last current
  // make new document the current (calls "slotChangeView()" indirectly)
  //WorkView->setCurrentTab(WorkView->tabAt(cNo));
  // Added by Shivam
   WorkView->setCurrentIndex(cNo); 
  // must be done here, in order to call the change slot !

  // if only an untitled document was open -> close it
  if(view->Docs.count() == 2)
    if(view->Docs.getFirst()->DocName.isEmpty())
      if(!view->Docs.getFirst()->DocChanged)
        view->Docs.removeRef(view->Docs.getFirst());

  d->sizeOfAll(d->UsedX1, d->UsedY1, d->UsedX2, d->UsedY2);
  view->viewport()->repaint();
  view->drawn = false;
  return true;
}

// #######################################################################
// Changes to the next document in the TabBar.
void QucsApp::slotNextTab()
{
  int No = view->Docs.at() + 1;
  if(view->Docs.current() == view->Docs.getLast())
    No = 0;

  // make new document the current (calls "slotChangeView()" indirectly)
  //WorkView->setCurrentTab(WorkView->tabAt(No));
    // Added by Shivam
         WorkView->setCurrentIndex(No);
  view->viewport()->repaint();
  view->drawn = false;
}

// #######################################################################
void QucsApp::slotFileSettings()
{
  view->editText->setHidden(true); // disable text edit of component property

  SettingsDialog *d = new SettingsDialog(view->Docs.current(), this);
  d->exec();
  view->viewport()->repaint();
  view->drawn = false;
}

// --------------------------------------------------------------
void QucsApp::slotApplSettings()
{
  view->editText->setHidden(true); // disable text edit of component property

  QucsSettingsDialog *d = new QucsSettingsDialog(this);
  d->exec();
  view->viewport()->repaint();
  view->drawn = false;
}

// --------------------------------------------------------------
void QucsApp::slotFileNew()
{
  statusBar()->message(tr("Creating new schematic..."));
  view->editText->setHidden(true); // disable text edit of component property

  view->Docs.append(new QucsDoc(this, ""));
  // make new document the current
  //WorkView->setCurrentTab(WorkView->tabAt(view->Docs.at()));
    // Added by Shivam
          WorkView->setCurrentIndex(view->Docs.at());

  statusBar()->message(tr("Ready."));
}

// --------------------------------------------------------------
void QucsApp::slotFileOpen()
{
  static QString lastDir;  // to remember last directory and file
  view->editText->setHidden(true); // disable text edit of component property

  statusBar()->message(tr("Opening file..."));

  QString s = Q3FileDialog::getOpenFileName(
	lastDir.isEmpty() ? QString(".") : lastDir,
	QucsFileFilter, this, "", tr("Enter a Schematic Name"));

  if(s.isEmpty())
    statusBar()->message(tr("Opening aborted"), 2000);
  else {
    gotoPage(s);
    lastDir = s;   // remember last directory and file
    statusBar()->message(tr("Ready."));
  }
}


// ########################################################################
void QucsApp::updatePortNumber(int No)
{
  if(No<0) return;
  QString pathName = view->Docs.current()->DocName;
  QFileInfo Info(pathName);
  QString tmpStr, File, Name = Info.fileName();

  // enter new port number into ListView
  Q3ListViewItem *p;
  for(p = ConSchematics->firstChild(); p!=0; p = p->nextSibling()) {
    if(p->text(0) == Name) {
      if(No == 0) p->setText(1,"");
      else p->setText(1,QString::number(No)+tr("-port"));
      break;
    }
  }
  if(No == 0) return;

  // update all occurencies of subcircuit in all open documents
  QucsDoc *d;
  int DocNo = view->Docs.at();
  for(d = view->Docs.first(); d!=0; d = view->Docs.next())
    for(Component *pc=d->Comps->first(); pc!=0; pc=d->Comps->next())
      if(pc->Model == "Sub") {
	File = pc->Props.getFirst()->Value;
	if((File == pathName) || (File == Name)) {
	  d->Comps->setAutoDelete(false);
	  d->deleteComp(pc);
	  tmpStr = pc->Name;
	  //((Subcircuit*)pc)->recreate();
	  pc->Name = tmpStr;
	  d->insertRawComponent(pc);
	  d->Comps->setAutoDelete(true);
	}
      }

  view->Docs.at(DocNo);  // back to the last current document
}


// ######################################################################
bool QucsApp::saveCurrentFile()
{
  if(view->Docs.current()->DocName.isEmpty())
    return saveAs();

  view->Docs.current()->PosX = view->contentsX();
  view->Docs.current()->PosY = view->contentsY();

  int Result = view->Docs.current()->save();   // SAVE
  if(Result < 0)  return false;

  updatePortNumber(Result);
  return true;
}

// ######################################################################
void QucsApp::slotFileSave()
{
  statusBar()->message(tr("Saving file..."));
  view->blockSignals(true);   // no user interaction during that time
  view->editText->setHidden(true); // disable text edit of component property

  if(!saveCurrentFile()) {
    view->blockSignals(false);
    statusBar()->message(tr("Saving aborted"), 2000);
    statusBar()->message(tr("Ready."));
    return;
  }

  view->blockSignals(false);
  statusBar()->message(tr("Ready."));
}

// #######################################################################
bool QucsApp::saveAs()
{
  static QString lastDir;  // to remember last directory and file

  int n = -1;
  bool intoView = true;
  QString s, Filter;
  QFileInfo Info;
  while(true) {
    intoView = true;
    s = view->Docs.current()->DocName;
    Info.setFile(s);
    if(s.isEmpty()) {   // which is default directory ?
      if(ProjName.isEmpty()) {
        if(lastDir.isEmpty())  s = QDir::currentDirPath();
        else  s = lastDir;
      }
      else s = QucsWorkDir.path();
    }

    s = Q3FileDialog::getSaveFileName(s,
	Info.extension() == "dpl" ? tr("Data Display")+" (*.dpl)"
				  : QucsFileFilter,
	this, "", tr("Enter a Document Name"), &Filter);
    if(s.isEmpty())  
      return false;
    Info.setFile(s);                 // try to guess the best extension ...
    if(Info.extension(false).isEmpty())  // ... if no one was specified
    {
      if(Filter.right(7) == "(*.dpl)")  
        s += ".dpl";
      else 
        s += ".sch";
    }

    Info.setFile(s);
    if(QFile::exists(s)) {
      n = QMessageBox::warning(this, tr("Warning"),
		tr("The file '")+Info.fileName()+tr("' already exists!\n")+
		tr("Saving will overwrite the old one! Continue?"),
		tr("No"), tr("Yes"), tr("Cancel"));
      if(n == 2) return false;    // cancel
      if(n == 0) continue;
      intoView = false;    // file already exists
    }

    QString ext = Info.extension(false);
    if(ext != "sch") if(ext != "dpl") {
      n = QMessageBox::warning(this, tr("Warning"),
		tr("Only the extensions '.sch' and '.dpl'\n")+
		tr("will appear in the content browser! Continue?"),
		tr("No"), tr("Yes"), tr("Cancel"));
      if(n == 2) return false;    // cancel
      if(n == 0) continue;
    }

    // search, if document is open
    QucsDoc *d;
    int No = view->Docs.at();
    for(d = view->Docs.first(); d!=0; d = view->Docs.next())
      if(d->DocName == s) break;
    view->Docs.at(No);   // back to old current document
    if(d) {
      QMessageBox::information(this, tr("Info"),
		tr("Cannot overwrite an open document"));
      return false;
    }

    break;
  }
  view->Docs.current()->setName(s);
  lastDir = Info.dirPath(true);  // remember last directory and file

  if(intoView) {    // insert new name in Content ListView ?
//    Info.setFile(s);
    if(Info.dirPath(true) == QucsWorkDir.absPath())
      if(!ProjName.isEmpty()) {
	s = Info.fileName();  // remove path from file name
	if(Info.extension(false) == "dpl")
	  Content->setSelected(new Q3ListViewItem(ConDisplays, s), true);
	else
	  Content->setSelected(new Q3ListViewItem(ConSchematics, s), true);
      }
  }

  view->Docs.current()->PosX = view->contentsX();
  view->Docs.current()->PosY = view->contentsY();
  n = view->Docs.current()->save();   // SAVE
  if(n < 0)  return false;

  updatePortNumber(n);
  return true;
}

// #######################################################################
void QucsApp::slotFileSaveAs()
{
  statusBar()->message(tr("Saving file under new filename..."));
  view->blockSignals(true);   // no user interaction during the time
  view->editText->setHidden(true); // disable text edit of component property

  if(!saveAs()) {
    view->blockSignals(false);
    statusBar()->message(tr("Saving aborted"), 3000);
    statusBar()->message(tr("Ready."));
    return;
  }

  view->blockSignals(false);
  statusBar()->message(tr("Ready."));
}


// #######################################################################
void QucsApp::slotFileSaveAll()
{
  statusBar()->message(tr("Saving all files..."));
  view->editText->setHidden(true); // disable text edit of component property

  QucsDoc *tmp = view->Docs.current();  // remember the current
  view->blockSignals(true);   // no user interaction during the time

  for(QucsDoc *ptr = view->Docs.first(); ptr != 0; ptr = view->Docs.next()) {
    if(ptr->DocName.isEmpty())  // make document the current ?
    //  WorkView->setCurrentTab(WorkView->tabAt(view->Docs.at()));
      // Added by Shivam
            WorkView->setCurrentIndex(view->Docs.at());
    saveCurrentFile();
  }

  view->Docs.findRef(tmp);
  view->blockSignals(false);
  // back to the current document
  //WorkView->setCurrentTab(WorkView->tabAt(view->Docs.at()));
    // Added by Shivam
          WorkView->setCurrentIndex(view->Docs.at());

  statusBar()->message(tr("Ready."));
}

// #######################################################################
void QucsApp::slotFileClose()
{
  statusBar()->message(tr("Closing file..."));
  view->editText->setHidden(true); // disable text edit of component property

  if(view->Docs.current()->DocChanged) {
    switch(QMessageBox::warning(this,tr("Closing ifREEDA document"),
      tr("The document contains unsaved changes!\n")+
      tr("Do you want to save the changes before closing?"),
      tr("&Save"), tr("&Discard"), tr("Cancel"), 0, 2)) {
      case 0 : slotFileSave();
               break;
      case 2 : return;
    }
  }

  view->Docs.remove();

  if(view->Docs.isEmpty())  // if no document left, create an untitled
    view->Docs.append(new QucsDoc(this, ""));

  nextDocument(true);

  QucsDoc *d = view->Docs.current();
  // which mode: schematic or symbol editor ?
  if((CompChoose->count() > 1) == d->symbolMode)
    changeSchematicSymbolMode(d);

  statusBar()->message(tr("Ready."));
}

// ######################################################################
void QucsApp::slotFilePrint()
{
  statusBar()->message(tr("Printing..."));
  view->editText->setHidden(true); // disable text edit of component property

  if (Printer->setup(this))  // print dialog
  {
    QPainter painter;
    painter.begin(Printer);
    view->Docs.current()->print(&painter, true);
    painter.end();
  };

  statusBar()->message(tr("Ready."));
}

// ######################################################################
void QucsApp::slotFilePrintSelected()
{
  statusBar()->message(tr("Printing..."));
  view->editText->setHidden(true); // disable text edit of component property

  if (Printer->setup(this))  // print dialog
  {
    QPainter painter;
    painter.begin(Printer);
    view->Docs.current()->print(&painter, false);
    painter.end();
  };

  statusBar()->message(tr("Ready."));
}

// --------------------------------------------------------------------
// Exits the application.
void QucsApp::slotFileQuit()
{
  statusBar()->message(tr("Exiting application..."));
  view->editText->setHidden(true); // disable text edit of component property

  int exit = QMessageBox::information(this, tr("Quit..."),
				      tr("Do you really want to quit?"),
				      tr("Yes"), tr("No"));
  if(exit == 0)
    if(closeAllFiles()) {
      emit signalKillEmAll();   // kill all subprocesses
      qApp->quit();
    }

  statusBar()->message(tr("Ready."));
}

//-----------------------------------------------------------------
// To get all close events.
void QucsApp::closeEvent(QCloseEvent* Event)
{
  Event->ignore();
  if(closeAllFiles()) {
    emit signalKillEmAll();   // kill all subprocesses
    qApp->quit();
  }
}

// --------------------------------------------------------------------
void QucsApp::slotEditCut()
{
  statusBar()->message(tr("Cutting selection..."));
  view->editText->setHidden(true); // disable text edit of component property

  QClipboard *cb = QApplication::clipboard();   // get system clipboard
  QString s = view->Docs.current()->copySelected(true);
  if(!s.isEmpty()) {
    cb->setText(s, QClipboard::Clipboard);
    view->viewport()->repaint();
  }

  statusBar()->message(tr("Ready."));
}

// ######################################################################
void QucsApp::slotEditCopy()
{
  statusBar()->message(tr("Copying selection to clipboard..."));
  view->editText->setHidden(true); // disable text edit of component property

  QClipboard *cb = QApplication::clipboard();   // get system clipboard
  QString s = view->Docs.current()->copySelected(false);
  if(!s.isEmpty())
    cb->setText(s, QClipboard::Clipboard);

  statusBar()->message(tr("Ready."));
}

// ########################################################################
// Is called when the toolbar button is pressed to go into a subcircuit.
void QucsApp::slotIntoHierarchy()
{
  view->editText->setHidden(true); // disable text edit of component property

  QucsDoc *Doc = view->Docs.current();
  Component *pc = view->Docs.current()->searchSelSubcircuit();
  if(pc == 0) return;

  QString *ps = new QString("*");
  HierarchyHistory.append(ps);    // sign to not clear HierarchyHistory

  QString s = QucsWorkDir.filePath(pc->Props.getFirst()->Value);
  if(!gotoPage(s)) return;

  *(HierarchyHistory.getLast()) = Doc->DocName; // remember for the way back
  popH->setEnabled(true);
}

// ########################################################################
// Is called when the toolbar button is pressed to leave a subcircuit.
void QucsApp::slotPopHierarchy()
{
  view->editText->setHidden(true); // disable text edit of component property

  if(HierarchyHistory.count() == 0) return;

  QString Doc = *(HierarchyHistory.getLast());
  *(HierarchyHistory.last()) = "*";  // sign to not clear HierarchyHistory

  if(!gotoPage(Doc)) return;

  HierarchyHistory.remove();
  if(HierarchyHistory.count() == 0) popH->setEnabled(false);
}

// ########################################################################
void QucsApp::slotShowAll()
{
  view->editText->setHidden(true); // disable text edit of component property

  int x1 = view->Docs.current()->UsedX1;
  int y1 = view->Docs.current()->UsedY1;
  int x2 = view->Docs.current()->UsedX2;
  int y2 = view->Docs.current()->UsedY2;

  if(x2==0) if(y2==0) if(x1==0) if(y1==0) return;
  x1 -= 40;  y1 -= 40;
  x2 += 40;  y2 += 40;

  float xScale = float(view->visibleWidth()) / float(x2-x1);
  float yScale = float(view->visibleHeight()) / float(y2-y1);
  if(xScale > yScale) xScale = yScale;
  if(xScale > 10.0) xScale = 10.0f;
  if(xScale < 0.01) xScale = 0.01f;
  view->Docs.current()->Scale = xScale;

  view->Docs.current()->ViewX1 = x1;
  view->Docs.current()->ViewY1 = y1;
  view->Docs.current()->ViewX2 = x2;
  view->Docs.current()->ViewY2 = y2;
  view->resizeContents(int(xScale*float(x2-x1)), int(xScale*float(y2-y1)));

  view->viewport()->repaint();
  view->drawn = false;
}

// -----------------------------------------------------------
// Sets the scale factor to 1.
void QucsApp::slotShowOne()
{
  view->editText->setHidden(true); // disable text edit of component property
  QucsDoc *d = view->Docs.current();

  d->Scale = 1.0;

  int x1 = d->UsedX1;
  int y1 = d->UsedY1;
  int x2 = d->UsedX2;
  int y2 = d->UsedY2;

//  view->Docs.current()->sizeOfAll(x1, y1, x2, y2);
  if(x2==0) if(y2==0) if(x1==0) if(y1==0) x2 = y2 = 800;

  d->ViewX1 = x1-40;
  d->ViewY1 = y1-40;
  d->ViewX2 = x2+40;
  d->ViewY2 = y2+40;
  view->resizeContents(x2-x1+80, y2-y1+80);
  view->viewport()->repaint();
  view->drawn = false;
}

void QucsApp::slotZoomOut()
{
  view->editText->setHidden(true); // disable text edit of component property
  
  view->Zoom(0.5);
  view->viewport()->repaint();
  view->drawn = false;
}

// -----------------------------------------------------------------------
// Is called when the simulate toolbar button is pressed.
void QucsApp::slotGenNetlist()
{
  view->Docs.current()->GenNetlist = true;//Only create a netlist, no simulation
  
  view->editText->setHidden(true); // disable text edit of component property
  
  if(view->Docs.current()->DocName.isEmpty()) // if document 'untitled' ...
    if(!saveCurrentFile()) return;            // ... save schematic before

  SimMessage *sim = new SimMessage(view->Docs.current(), this);
  // disconnect is automatically performed, if one of the involved objects
  // is destroyed !
  connect(sim, SIGNAL(SimulationEnded(int, SimMessage*)), this,
		SLOT(slotAfterSimulation(int, SimMessage*)));
  connect(sim, SIGNAL(displayDataPage(QString)),
		this, SLOT(slotChangePage(QString)));

  sim->show();
  if(!sim->startProcess()) return;

  // to kill it before qucs ends
  connect(this, SIGNAL(signalKillEmAll()), sim, SLOT(slotClose()));
}

// -----------------------------------------------------------------------
// Is called when the simulate toolbar button is pressed.
void QucsApp::slotSimulate()
{
  view->Docs.current()->GenNetlist = false;//Only create a netlist, no simulation

  view->editText->setHidden(true); // disable text edit of component property
  
  if(view->Docs.current()->DocName.isEmpty()) // if document 'untitled' ...
    if(!saveCurrentFile()) return;            // ... save schematic before

  SimMessage *sim = new SimMessage(view->Docs.current(), this);
  // disconnect is automatically performed, if one of the involved objects
  // is destroyed !
  connect(sim, SIGNAL(SimulationEnded(int, SimMessage*)), this,
		SLOT(slotAfterSimulation(int, SimMessage*)));
  connect(sim, SIGNAL(displayDataPage(QString)),
		this, SLOT(slotChangePage(QString)));

  sim->show();
  if(!sim->startProcess()) return;

  // to kill it before qucs ends
  connect(this, SIGNAL(signalKillEmAll()), sim, SLOT(slotClose()));
}

// ------------------------------------------------------------------------
// Is called after the simulation process terminates.
void QucsApp::slotAfterSimulation(int Status, SimMessage *sim)
{
  bool shouldClosed = false; // simulation window close after simualtion ?

  if(Status == 0) {  // no errors ocurred ?
    if(view->Docs.current()->SimOpenDpl) {
      shouldClosed = true;
      //slotChangePage(sim->Doc->DataDisplay);  // switch to data display
      view->viewport()->update();
    }
    else {
      view->Docs.current()->reloadGraphs();  // load recent simulation data
      view->viewport()->update();
    }

    // put all dataset files into "Content"-ListView (update)
/*    QStringList Elements = ProjDir.entryList("*.dat", QDir::Files, QDir::Name);
    for(it = Elements.begin(); it != Elements.end(); ++it)
      new QListViewItem(ConDatasets, (*it).ascii());*/
  }

  if(shouldClosed) sim->slotClose();  // close and delete simulation window
}


/***************************************************
*
*            Possible code to convert freeda data file
*            into format ifreeda can display need to 
*            determine qucs format of data file
*
********************************************************
  QFile file;
  QString fileName = "netlist.net";
  file.setName(fileName);

  file.setName(QucsWorkDir.filePath(file.name()));
  if(!file.open(IO_ReadOnly)) {
    QMessageBox::critical(0, QObject::tr("Error"),
                 QObject::tr("Cannot find model file: ")+file.name());
  }
  // *****************************************************************
  // To strongly speed up the file read operation the whole file is
  // read into the memory in one piece.
  QTextStream ReadWhole(&file);
  QString FileString = ReadWhole.read();
  file.close();

  QString Line;

  int i=0, j=0, length;   // i and j must not be used temporarily !!!

  i = FileString.find(".out");//Find the start of the .out statements.
  Line = FileString.remove(0, i);
  j = Line.find(".end");//Find the end of the .out parameters.
  Line.truncate(j);

  Line = Line.simplifyWhiteSpace();
  QString path = Line.section(' ',-1);
  QString dataFile = path.section('/',-1);
  length = dataFile.length();
  dataFile = dataFile.remove(length-1, 1);
  printf(path);
  printf("\n");
  printf(dataFile);
  
  QFile freedaDataFile;
  fileName = dataFile;
  freedaDataFile.setName(fileName);

  freedaDataFile.setName(QucsWorkDir.filePath(freedaDataFile.name()));
  if(!freedaDataFile.open(IO_ReadOnly)) {
    QMessageBox::critical(0, QObject::tr("Error"),
                 QObject::tr("Cannot find model file: ")+freedaDataFile.name());
  }
  // *****************************************************************
  // To strongly speed up the file read operation the whole file is
  // read into the memory in one piece.
  QTextStream ReadWholeFreeda(&freedaDataFile);
  FileString = ReadWholeFreeda.read();
  freedaDataFile.close();
 
//  QVector<int> x;
//  QVector<int> y;

    QString x, y;

  x = "<fREEDA Dataset>\n<indep x current>\n";
  y = "<dep y output>\n";     
  FileString = FileString.simplifyWhiteSpace();
  length = FileString.length();
  for(i=0; i<=length; i++)
  {
     if (i%2 ==0)
       x = x + FileString.section(' ', i, i) + "\n";
     else 
       y = y + FileString.section(' ', i, i) + "\n";
  }
  x = x.simplifyWhiteSpace();
  y = y.simplifyWhiteSpace();  
  
  x = x + "</indep>\n";
  y = y + "</dep>\n";

  QFile         ifreedaDataFile;
  QTextStream   stream;
  
  ifreedaDataFile.setName(QucsWorkDir.filePath(dataFile+".dat"));
  if(!ifreedaDataFile.open(IO_WriteOnly)) {
    //ErrText->insert(tr("ERROR: Cannot write ifREEDA Data File file!"));
  }
  stream.setDevice(&ifreedaDataFile);
  stream << x << y;
  ifreedaDataFile.close();
***********************************************/
// ------------------------------------------------------------------------
// Changes to the corresponding data display page or vice versa.
void QucsApp::slotChangePage(QString Name)
{
  if(Name.isEmpty())  return;

  QucsDoc   *d;
  QFileInfo Info;
  QucsDoc   *Doc = view->Docs.current();

  // search, if page is already loaded
  for(d = view->Docs.first(); d!=0; d = view->Docs.next()) {
    Info.setFile(d->DocName);
    if(Info.fileName() == Name) break;
  }

  if(d == 0) {   // no open page found ?
    Info.setFile(Name);
    if(Info.isRelative())
      Name = QucsWorkDir.filePath(Name);
    QFile file(Name);
    if(!file.open(QIODevice::ReadOnly))       // load page
    {
      if(!file.open(QIODevice::ReadWrite)) 
      {  // if it doesn't exist, create
        view->Docs.findRef(Doc);
        QMessageBox::critical(this, tr("Error"), tr("Cannot create ")+Info.dirPath(true)+
    			QDir::convertSeparators ("/")+Name);
        return;
      }
      else 
        new Q3ListViewItem(ConDisplays, Info.fileName()); // add new name
    }
    file.close();

    d = new QucsDoc(this, Name);
    view->Docs.append(d);   // create new page

    if(!d->load()) 
    {
      view->Docs.remove();
      view->Docs.findRef(Doc);
      view->viewport()->repaint();
      view->drawn = false;
      return;
    }
  }

  int cNo = view->Docs.at();
  view->Docs.findRef(Doc);  // back to the last current
  //if(WorkView->indexOf(WorkView->currentTab()) == cNo)  // if page not ...
  // Added by Shivam In Qt4.3+ id and index means same.
    if(WorkView->currentIndex() == cNo)
    view->Docs.current()->reloadGraphs();       // ... changes, reload here !
  //WorkView->setCurrentTab(WorkView->tabAt(cNo));  // make new doc the current
  //Added by Shivam
   WorkView->setCurrentIndex(cNo);   
  // must be done before the next step, in order to call the change slot !

  d->sizeOfAll(d->UsedX1, d->UsedY1, d->UsedX2, d->UsedY2);

  TabView->setCurrentPage(2);   // switch to "Component"-Tab
  if(Name.right(4) == ".dpl") {
    CompChoose->setCurrentItem(COMBO_Diagrams);   // switch to diagrams
    slotSetCompView(COMBO_Diagrams);
  }
}

// ------------------------------------------------------------------------
// Changes to the data display of current page.
void QucsApp::slotToPage()
{
  QString Name = view->Docs.current()->DataDisplay;
  if(Name.isEmpty()) {
    QMessageBox::critical(this, tr("Error"), tr("No page set !"));
    return;
  }

  slotChangePage(Name);
}

// #######################################################################
// Is called when a double click is made in the content ListView.
//void QucsApp::slotOpenContent(QListViewItem *item, const QPoint &, int column)   // QT 3.2
void QucsApp::slotOpenContent(Q3ListViewItem *item)
{
  view->editText->setHidden(true); // disable text edit of component property

  if(item == 0) return;   // no item was double clicked
  if(item->parent() == 0) return;  // no component, but item "schematic", ...

  QucsWorkDir.setPath(QucsHomeDir.path());
  QString p = ProjName+"_prj";
  if(!QucsWorkDir.cd(p)) {
    QMessageBox::critical(this, tr("Error"),
                          tr("Cannot access project directory: ")+
			  QucsHomeDir.path()+QDir::convertSeparators ("/")+p);
    return;
  }

  QFileInfo Info(QucsWorkDir.filePath(item->text(0)));
  if(Info.extension(false) == "dat") {
    Acts.editFile(Info.absFilePath());  // open datasets with text editor
    return;
  }
  gotoPage(Info.absFilePath());

  if(item->text(1).isEmpty()) return;
  // switch on the 'select' action
  Acts.select->blockSignals(true);
  Acts.select->setOn(true);
  Acts.select->blockSignals(false);
  Acts.slotSelect(true);
}

// ########################################################################
// Is called when the close project menu is called.
void QucsApp::slotMenuCloseProject()
{
  view->editText->setHidden(true); // disable text edit of component property

  if(!closeAllFiles()) return;   // close files and ask for saving them
  QucsDoc *d = new QucsDoc(this, "");
  view->Docs.append(d);   // create 'untitled' file
  view->resizeContents(int(d->Scale*double(d->ViewX2-d->ViewX1)),
                       int(d->Scale*double(d->ViewY2-d->ViewY1)));
  view->setContentsPos(d->PosX, d->PosY);
  view->viewport()->repaint();
  view->drawn = false;

  setCaption("fREEDA " PACKAGE_VERSION + tr(" - Project: "));
  QucsWorkDir.setPath(QDir::homeDirPath()+QDir::convertSeparators ("/freeda"));

  Content->setColumnText(0,tr("Content of")); // column text

  Content->clear();   // empty content view
  ConDatasets   = new Q3ListViewItem(Content, tr("Datasets"));
  ConDisplays   = new Q3ListViewItem(Content, tr("Data Displays"));
  ConSchematics = new Q3ListViewItem(Content, tr("Schematics"));

  TabView->setCurrentPage(0);   // switch to "Projects"-Tab
  ProjName = "";
}

// ########################################################################
// Is called when the open project menu is called.
void QucsApp::slotMenuOpenProject()
{
  Q3FileDialog *d = new Q3FileDialog(QucsHomeDir.path());
  d->setCaption(tr("Choose Project Directory for Opening"));
  d->setShowHiddenFiles(true);
  d->setMode(Q3FileDialog::DirectoryOnly);
  if(d->exec() != QDialog::Accepted) return;

  QString s = d->selectedFile();
  if(s.isEmpty()) return;

  s = s.left(s.length()-1);   // cut off trailing '/'
  int i = s.findRev('/');
  if(i > 0) s = s.mid(i+1);   // cut out the last subdirectory
  s.remove("_prj");
  OpenProject(d->selectedFile(), s);
}

// #######################################################################
// Is called when the open project button is pressed.
void QucsApp::slotOpenProject(Q3ListBoxItem *item)
{
  OpenProject(QucsHomeDir.filePath(item->text()+"_prj"), item->text());
}

// ########################################################################
// Checks whether this file is a qucs file and whether it is an subcircuit.
// It returns the number of subcircuit ports.
int QucsApp::testFile(const QString& DocName)
{
  QFile file(DocName);
  if(!file.open(QIODevice::ReadOnly)) {
    return -1;
  }

  QString Line;
  // *****************************************************************
  // To strongly speed up the file read operation the whole file is
  // read into the memory in one piece.
  Q3TextStream ReadWhole(&file);
  QString FileString = ReadWhole.read();
  file.close();
  Q3TextStream stream(&FileString, QIODevice::ReadOnly);


  // read header **************************
  do {
    if(stream.atEnd()) {
      file.close();
      return -2;
    }
    Line = stream.readLine();
    Line = Line.stripWhiteSpace();
  } while(Line.isEmpty());
  
  if(Line.left(16) != "<Qucs Schematic ") {  // wrong file type ?
    file.close();
    return -3;
  }

  QString s = PACKAGE_VERSION;
  s.remove('.');
  Line = Line.mid(16, Line.length()-17);
  Line.remove('.');
  if(Line > s) {  // wrong version number ? (only backward compatible)
    file.close();
    return -4;
  }

  // read content *************************
  while(!stream.atEnd()) {
    Line = stream.readLine();
    if(Line == "<Components>") break;
  }

  int z=0;
  while(!stream.atEnd()) {
    Line = stream.readLine();
    if(Line == "</Components>") {
      file.close();
      return z;       // return number of ports
    }

    Line = Line.stripWhiteSpace();
    s    = Line.section(' ',0,0);    // component type
    if(s == "<Port") z++;
  }
  return -5;  // component field not closed
}

// #######################################################################
// Opens an existing project.
void QucsApp::OpenProject(const QString& Path, const QString& Name)
{
  view->editText->setHidden(true); // disable text edit of component property

  if(!closeAllFiles()) return;   // close files and ask for saving them
  QucsDoc *d = new QucsDoc(this, "");
  view->Docs.append(d);   // create 'untitled' file
  view->resizeContents(int(d->Scale*double(d->ViewX2-d->ViewX1)),
                       int(d->Scale*double(d->ViewY2-d->ViewY1)));
  view->setContentsPos(d->PosX, d->PosY);
  view->viewport()->repaint();
  view->drawn = false;

  QDir ProjDir(QDir::cleanDirPath(Path));
  if(!ProjDir.exists() || !ProjDir.isReadable()) { // check project directory

    QMessageBox::critical(this, tr("Error"),
                          tr("Cannot access project directory : ")+Path);
    return;
  }
  QucsWorkDir.setPath(ProjDir.path());

  Content->setColumnText(0,tr("Content of '")+Name+tr("'")); // column text
//  Content->setColumnWidth(0, Content->width()-5);

  Content->clear();   // empty content view
  ConDatasets   = new Q3ListViewItem(Content, tr("Datasets"));
  ConDisplays   = new Q3ListViewItem(Content, tr("Data Displays"));
  ConSchematics = new Q3ListViewItem(Content, tr("Schematics"));

  int n;
  // put all schematic files into "Content"-ListView
  QStringList Elements = ProjDir.entryList("*.sch", QDir::Files, QDir::Name);
  QStringList::iterator it;
  for(it = Elements.begin(); it != Elements.end(); ++it) {
    n = testFile(ProjDir.filePath((*it).ascii()));
    if(n >= 0) {
      if(n > 0)
        new Q3ListViewItem(ConSchematics, (*it).ascii(),
			  QString::number(n)+tr("-port"));
      else new Q3ListViewItem(ConSchematics, (*it).ascii());
    }
  }

  // put all data display files into "Content"-ListView
  Elements = ProjDir.entryList("*.dpl", QDir::Files, QDir::Name);
  for(it = Elements.begin(); it != Elements.end(); ++it)
    new Q3ListViewItem(ConDisplays, (*it).ascii());

  // put all dataset files into "Content"-ListView
  Elements = ProjDir.entryList("*.dat", QDir::Files, QDir::Name);
  for(it = Elements.begin(); it != Elements.end(); ++it)
    new Q3ListViewItem(ConDatasets, (*it).ascii());

  TabView->setCurrentPage(1);   // switch to "Content"-Tab
  Content->firstChild()->setOpen(true);  // show schematics
  ProjName = Name;   // remember the name of project

  // show name in title of main window
  setCaption("ifREEDA 1.0" + tr(" - Project: ")+Name);
}

// #######################################################################
// Is called, when "Create New Project"-Button is pressed.
void QucsApp::slotProjNewButt()
{
  view->editText->setHidden(true); // disable text edit of component property

  NewProjDialog *d = new NewProjDialog(this);
  if(d->exec() != QDialog::Accepted) return;

  QDir projDir(QucsHomeDir.path());  
  if(projDir.mkdir(d->ProjName->text()+"_prj")) {
    Projects->insertItem(d->ProjName->text(),0);  // at first position
    if(d->OpenProj->isChecked()) slotOpenProject(Projects->firstItem());
  }
  else QMessageBox::information(this, tr("Info"),
                    tr("Cannot create project directory!"));
}

// ######################################################################
// The following arrays contains the components that appear in the
// component listview.
typedef Component*  (*pInfoFunc) (QString&, char* &, bool);
pInfoFunc Simulations[] =
  {0};

pInfoFunc lumpedComponents[] =
  {0};

pInfoFunc Sources[] =
  {0};

pInfoFunc TransmissionLines[] =
  {0};

pInfoFunc nonlinearComps[] =
  {0};

// #######################################################################
// Whenever the Component Library ComboBox is changed, this slot fills the
// Component IconView with the appropriat components.
void QucsApp::slotSetCompView(int index)
{
  view->editText->setHidden(true); // disable text edit of component property

  char *File;
  int i = 0;
  QString Name;
  pInfoFunc *Infos = 0;

  CompComps->clear();   // clear the IconView
  if((index+1) >= CompChoose->count()) {
    //new Q3IconViewItem(CompComps, tr("Line"),QImage(QucsSettings.BitmapDir + "line.png"));
    //new Q3IconViewItem(CompComps, tr("Arrow"),QImage(QucsSettings.BitmapDir + "arrow.png"));
    //new Q3IconViewItem(CompComps, tr("Text"),QImage(QucsSettings.BitmapDir + "text.png"));
    //new Q3IconViewItem(CompComps, tr("Ellipse"),QImage(QucsSettings.BitmapDir + "ellipse.png"));
    //new Q3IconViewItem(CompComps, tr("Rectangle"),QImage(QucsSettings.BitmapDir + "rectangle.png"));
    //new Q3IconViewItem(CompComps, tr("filled Ellipse"),QImage(QucsSettings.BitmapDir + "filledellipse.png"));
    //new Q3IconViewItem(CompComps, tr("filled Rectangle"),QImage(QucsSettings.BitmapDir + "filledrect.png"));
    //new Q3IconViewItem(CompComps, tr("Elliptic Arc"),QImage(QucsSettings.BitmapDir + "ellipsearc.png"));

 //Added by Shivam
    new Q3IconViewItem(CompComps, tr("Line"),QPixmap((QucsSettings.BitmapDir + "line.png"),"PNG",0));
    new Q3IconViewItem(CompComps, tr("Arrow"),QPixmap((QucsSettings.BitmapDir + "arrow.png"),"PNG",0));
    new Q3IconViewItem(CompComps, tr("Text"),QPixmap((QucsSettings.BitmapDir + "text.png"),"PNG",0));
    new Q3IconViewItem(CompComps, tr("Ellipse"),QPixmap((QucsSettings.BitmapDir + "ellipse.png"),"PNG",0));
    new Q3IconViewItem(CompComps, tr("Rectangle"),QPixmap((QucsSettings.BitmapDir + "rectangle.png"),"PNG",0));
    new Q3IconViewItem(CompComps, tr("filled Ellipse"),QPixmap((QucsSettings.BitmapDir + "filledellipse.png"),"PNG",0));
    new Q3IconViewItem(CompComps, tr("filled Rectangle"),QPixmap((QucsSettings.BitmapDir + "filledrect.png"),"PNG",0));
    new Q3IconViewItem(CompComps, tr("Elliptic Arc"),QPixmap((QucsSettings.BitmapDir + "ellipsearc.png"),"PNG",0));

    return;
  }

  switch(index) {
    //case COMBO_passive:   Infos = &lumpedComponents[0];  break;
    //case COMBO_Sources:   Infos = &Sources[0];           break;
    //case COMBO_TLines:    Infos = &TransmissionLines[0]; break;
    //case COMBO_nonlinear: Infos = &nonlinearComps[0];    break;
    case COMBO_File:
      //new Q3IconViewItem(CompComps, tr("SPICE netlist"),QImage(QucsSettings.BitmapDir + "spicefile.png"));
      //new Q3IconViewItem(CompComps, tr("1-port S parameter file"),QImage(QucsSettings.BitmapDir + "spfile1.png"));
      //new Q3IconViewItem(CompComps, tr("2-port S parameter file"),QImage(QucsSettings.BitmapDir + "spfile2.png"));
      //new Q3IconViewItem(CompComps, tr("3-port S parameter file"),QImage(QucsSettings.BitmapDir + "spfile3.png"));
      //new Q3IconViewItem(CompComps, tr("4-port S parameter file"),QImage(QucsSettings.BitmapDir + "spfile4.png"));
      //new Q3IconViewItem(CompComps, tr("5-port S parameter file"),QImage(QucsSettings.BitmapDir + "spfile5.png"));
      //new Q3IconViewItem(CompComps, tr("6-port S parameter file"),QImage(QucsSettings.BitmapDir + "spfile6.png"));

      //Added by Shivam
      new Q3IconViewItem(CompComps, tr("SPICE netlist"),QPixmap((QucsSettings.BitmapDir + "spicefile.png"),"PNG",0));
      new Q3IconViewItem(CompComps, tr("1-port S parameter file"),QPixmap((QucsSettings.BitmapDir + "spfile1.png"),"PNG",0));
      new Q3IconViewItem(CompComps, tr("2-port S parameter file"),QPixmap((QucsSettings.BitmapDir + "spfile2.png"),"PNG",0));
      new Q3IconViewItem(CompComps, tr("3-port S parameter file"),QPixmap((QucsSettings.BitmapDir + "spfile3.png"),"PNG",0));
      new Q3IconViewItem(CompComps, tr("4-port S parameter file"),QPixmap((QucsSettings.BitmapDir + "spfile4.png"),"PNG",0));
      new Q3IconViewItem(CompComps, tr("5-port S parameter file"),QPixmap((QucsSettings.BitmapDir + "spfile5.png"),"PNG",0));
      new Q3IconViewItem(CompComps, tr("6-port S parameter file"),QPixmap((QucsSettings.BitmapDir + "spfile6.png"),"PMG",0));

      return;
    //case COMBO_Sims:     Infos = &Simulations[0];  break;
    case COMBO_Diagrams:
      new Q3IconViewItem(CompComps, tr("Cartesian"),QPixmap((QucsSettings.BitmapDir + "rect.png"),"PNG",0));
      new Q3IconViewItem(CompComps, tr("Polar"),QPixmap((QucsSettings.BitmapDir + "polar.png"),"PNG",0));
      new Q3IconViewItem(CompComps, tr("Tabular"),QPixmap((QucsSettings.BitmapDir + "tabular.png"),"PNG",0));
      new Q3IconViewItem(CompComps, tr("Smith Chart"),QPixmap((QucsSettings.BitmapDir + "smith.png"),"PNG",0));
      new Q3IconViewItem(CompComps, tr("Admittance Smith"),QPixmap((QucsSettings.BitmapDir + "ysmith.png"),"PNG",0));
      new Q3IconViewItem(CompComps, tr("Polar-Smith Combi"),QPixmap((QucsSettings.BitmapDir + "polarsmith.png"),"PNG",0));
      new Q3IconViewItem(CompComps, tr("Smith-Polar Combi"),QPixmap((QucsSettings.BitmapDir + "smithpolar.png"),"PNG",0));
      new Q3IconViewItem(CompComps, tr("3D-Cartesian"),QPixmap((QucsSettings.BitmapDir + "rect3d.png"),"PNG",0));
      new Q3IconViewItem(CompComps, tr("Locus Curve"),QPixmap((QucsSettings.BitmapDir + "curve.png"),"PNG",0));
      return;

      case COMBO_fREEDA_lumped: //fREEDA-EDIT
      {
 #include "fREEDA_lumped.cpp"
        i = 0;
        while(strcmp(fREEDA_lumped[i], "0"))
	{ 
	  new Q3IconViewItem(CompComps, tr(fREEDA_lumped[i]),QPixmap((QucsSettings.BitmapDir+QString(fREEDA_lumped[i])+".png"),"PNG",0));
	  i++;
	}
        return;
      }
      
      case COMBO_fREEDA_nonlinear: //fREEDA-EDIT
      {
 #include "fREEDA_nonlinear.cpp"
        i = 0;
        while(strcmp(fREEDA_nonlinear[i], "0"))
	{ 
	  new Q3IconViewItem(CompComps, tr(fREEDA_nonlinear[i]),QPixmap((QucsSettings.BitmapDir+QString(fREEDA_nonlinear[i])+".png"),"PNG",0));
	  i++;
	}
        return;
      }
      
      case COMBO_fREEDA_sources: //fREEDA-EDIT
      {
 #include "fREEDA_sources.cpp"
        i = 0;
        while(strcmp(fREEDA_sources[i], "0"))
	{ 
	  new Q3IconViewItem(CompComps, tr(fREEDA_sources[i]),QPixmap((QucsSettings.BitmapDir+QString(fREEDA_sources[i])+".png"),"PNG",0));
	  i++;
	}
        return;
      }
      
      case COMBO_fREEDA_thermal: //fREEDA-EDIT
      {
 #include "fREEDA_thermal.cpp"
        int i = 0;
        while(strcmp(fREEDA_thermal[i], "0"))
	{ 
	  new Q3IconViewItem(CompComps, tr(fREEDA_thermal[i]), QPixmap((QucsSettings.BitmapDir+QString(fREEDA_thermal[i])+".png"),"PNG",0));
	  i++;
	}
        return;
      }
      
      case COMBO_fREEDA_logic:
      {
 #include "fREEDA_logic.cpp"
        i = 0;
        while(strcmp(fREEDA_logic[i], "0"))
	{ 
	  new Q3IconViewItem(CompComps, tr(fREEDA_logic[i]),QPixmap((QucsSettings.BitmapDir+QString(fREEDA_logic[i])+".png"),"PNG",0));
	  i++;
	}
        return;
      }  
             
      case COMBO_fREEDA_other: //fREEDA-EDIT
      {
 #include "fREEDA_other.cpp"
        i = 0;
        while(strcmp(fREEDA_other[i], "0"))
	{ 
	  new Q3IconViewItem(CompComps, tr(fREEDA_other[i]),QPixmap((QucsSettings.BitmapDir+QString(fREEDA_other[i])+".png"),"PNG",0));
	  i++;
	}
        return;
      }

      case COMBO_fanalysis: //fREEDA-EDIT
      {
 #include "fREEDA_analysis.cpp"
        i = 0;
        while(strcmp(fREEDA_analysis[i], "0"))
	{ 
	  new Q3IconViewItem(CompComps, tr(fREEDA_analysis[i]),QPixmap((QucsSettings.BitmapDir+QString(fREEDA_analysis[i])+".png"),"PNG",0));
	  i++;
	}
        return;
      }
      default:
         return;
  }

  while(*Infos != 0) {
    (**Infos) (Name, File, false);
    new Q3IconViewItem(CompComps, Name,QPixmap((QucsSettings.BitmapDir+QString(File)+".png"),"PNG",0));
    Infos++;
  }
}

// ----------------------------------------------------------------------
// Is called when the mouse is clicked within the Component QIconView.
void QucsApp::slotSelectComponent(Q3IconViewItem *item)
{
  view->editText->setHidden(true); // disable text edit of component property

  // delete previously selected elements
  if(view->selComp != 0)  delete view->selComp;
  if(view->selDiag != 0)  delete view->selDiag;
  if(view->selPaint != 0) delete view->selDiag;
  view->selComp  = 0;   // no component selected
  view->selDiag  = 0;   // no diagram selected
  view->selPaint = 0;   // no painting selected

  if(item == 0) {   // mouse button pressed not over an item ?
//    CompComps->setSelected(CompComps->currentItem(), false);  // deselect component in ViewList
    CompComps->clearSelection();  // deselect component in ViewList
    if(view->drawn) view->viewport()->repaint();
    view->drawn = false;
    return;
  }

  // toggle last toolbar button off
  if(activeAction) {
    activeAction->blockSignals(true); // do not call toggle slot
    activeAction->setOn(false);       // set last toolbar button off
    activeAction->blockSignals(false);
  }
  activeAction = 0;


  if((CompChoose->currentItem()+1) >= CompChoose->count()) {
    switch(CompComps->index(item)) { // in normal ("edit schematic") mode
	case 0: view->selPaint = new GraphicLine();   break;
	case 1: view->selPaint = new Arrow();         break;
	case 2: view->selPaint = new GraphicText();   break;
	case 3: view->selPaint = new class Ellipse();       break;
	case 4: view->selPaint = new class Rectangle();     break;
	case 5: view->selPaint = new class Ellipse(true);   break;
	case 6: view->selPaint = new class Rectangle(true); break;
	case 7: view->selPaint = new EllipseArc();    break;
    }

    if(view->drawn) view->viewport()->repaint();
    view->drawn = false;
    view->MouseMoveAction = &QucsView::MMovePainting;
    view->MousePressAction = &QucsView::MPressPainting;
    view->MouseReleaseAction = &QucsView::MouseDoNothing;
    view->MouseDoubleClickAction = &QucsView::MouseDoNothing;
    return;
  }

  pInfoFunc Infos = 0;
  switch(CompChoose->currentItem()) {
   /*case COMBO_passive:
	 Infos = lumpedComponents[CompComps->index(item)];
	 break;
    case COMBO_Sources:
	 Infos = Sources[CompComps->index(item)];
	 break;
    case COMBO_TLines:
	 Infos = TransmissionLines[CompComps->index(item)];
	 break;
    case COMBO_nonlinear:
	 Infos = nonlinearComps[CompComps->index(item)];
	 break;*/
    case COMBO_File:
         //if(CompComps->index(item) == 0)
	   //view->selComp = new SpiceFile();
	 //else
	   //view->selComp = new SParamFile(CompComps->index(item));
	 break;
    //case COMBO_Sims:
    //	 Infos = Simulations[CompComps->index(item)];
    //	 break;
    case COMBO_Diagrams:
          switch(CompComps->index(item)) {
              case 0: view->selDiag = new RectDiagram();  break;
              case 1: view->selDiag = new PolarDiagram(); break;
              case 2: view->selDiag = new TabDiagram();   break;
              case 3: view->selDiag = new SmithDiagram(); break;
              case 4: view->selDiag = new SmithDiagram(0,0,false); break;
              case 5: view->selDiag = new PSDiagram();  break;
              case 6: view->selDiag = new PSDiagram(0,0,false); break;
              case 7: view->selDiag = new Rect3DDiagram(); break;
              case 8: view->selDiag = new CurveDiagram();  break;
          }

          if(view->drawn) view->viewport()->repaint();
          view->drawn = false;
          view->MouseMoveAction = &QucsView::MMoveDiagram;
          view->MousePressAction = &QucsView::MPressDiagram;
          view->MouseReleaseAction = &QucsView::MouseDoNothing;
          view->MouseDoubleClickAction = &QucsView::MouseDoNothing;
          return;
	  
     case COMBO_fREEDA_lumped:
     {
#include "fREEDA_lumped.cpp"
     QString elem_type = fREEDA_lumped[CompComps->index(item)];
     //const char* name = elem_type.data();
     view->selComp = new Interface(elem_type);
     if(view->drawn) view->viewport()->repaint();
     view->drawn = false;
     view->MouseMoveAction = &QucsView::MMoveComponent;
     view->MousePressAction = &QucsView::MPressComponent;
     view->MouseReleaseAction = &QucsView::MouseDoNothing;
     view->MouseDoubleClickAction = &QucsView::MouseDoNothing;
     return;
     }

     case COMBO_fREEDA_nonlinear:
     {
#include "fREEDA_nonlinear.cpp"
     QString elem_type = fREEDA_nonlinear[CompComps->index(item)];
     //const char* name = elem_type.data();
     view->selComp = new Interface(elem_type);
     if(view->drawn) view->viewport()->repaint();
     view->drawn = false;
     view->MouseMoveAction = &QucsView::MMoveComponent;
     view->MousePressAction = &QucsView::MPressComponent;
     view->MouseReleaseAction = &QucsView::MouseDoNothing;
     view->MouseDoubleClickAction = &QucsView::MouseDoNothing;
     return;
     }
     
     case COMBO_fREEDA_sources:
     {
#include "fREEDA_sources.cpp"
     QString elem_type = fREEDA_sources[CompComps->index(item)];
     //const char* name = elem_type.data();
     view->selComp = new Interface(elem_type);
     if(view->drawn) view->viewport()->repaint();
     view->drawn = false;
     view->MouseMoveAction = &QucsView::MMoveComponent;
     view->MousePressAction = &QucsView::MPressComponent;
     view->MouseReleaseAction = &QucsView::MouseDoNothing;
     view->MouseDoubleClickAction = &QucsView::MouseDoNothing;
     return;
     }
     
     case COMBO_fREEDA_thermal:
     {
#include "fREEDA_thermal.cpp"
     QString elem_type = fREEDA_thermal[CompComps->index(item)];
     //const char* name = elem_type.data();
     view->selComp = new Interface(elem_type);
     if(view->drawn) view->viewport()->repaint();
     view->drawn = false;
     view->MouseMoveAction = &QucsView::MMoveComponent;
     view->MousePressAction = &QucsView::MPressComponent;
     view->MouseReleaseAction = &QucsView::MouseDoNothing;
     view->MouseDoubleClickAction = &QucsView::MouseDoNothing;
     return;
     }
     
     case COMBO_fREEDA_logic:
     {
#include "fREEDA_logic.cpp"
     QString elem_type = fREEDA_logic[CompComps->index(item)];
     //const char* name = elem_type.data();
     view->selComp = new Interface(elem_type);
     if(view->drawn) view->viewport()->repaint();
     view->drawn = false;
     view->MouseMoveAction = &QucsView::MMoveComponent;
     view->MousePressAction = &QucsView::MPressComponent;
     view->MouseReleaseAction = &QucsView::MouseDoNothing;
     view->MouseDoubleClickAction = &QucsView::MouseDoNothing;
     return;
     }
        
     case COMBO_fREEDA_other:
     {
#include "fREEDA_other.cpp"
     QString elem_type = fREEDA_other[CompComps->index(item)];
     //string type = elem_type.data();  Commented by Shivam as .data returns to Pointer to data stored in QString. Qchar* data() 

     if(elem_type == "ground")  //(type == "ground")  Commented by Shivam
     {
     	view->selComp = new Ground();
     }
     else if(elem_type == "localreference")   // (type == "localreference")  Commented by Shivam
     {
        view->selComp = new LocalReference();
     }
     else
     {
        view->selComp = new Interface(elem_type);
     }
     if(view->drawn) view->viewport()->repaint();
     view->drawn = false;
     view->MouseMoveAction = &QucsView::MMoveComponent;
     view->MousePressAction = &QucsView::MPressComponent;
     view->MouseReleaseAction = &QucsView::MouseDoNothing;
     view->MouseDoubleClickAction = &QucsView::MouseDoNothing;
     return;
     }
  
     case COMBO_fanalysis:
     {  
#include "fREEDA_analysis.cpp"
     QString analysis_type = fREEDA_analysis[CompComps->index(item)];
     //const char* name = analysis_type.data();
     view->selComp = new Interface(analysis_type);
     if(view->drawn) view->viewport()->repaint();
     view->drawn = false;
     view->MouseMoveAction = &QucsView::MMoveComponent;
     view->MousePressAction = &QucsView::MPressComponent;
     view->MouseReleaseAction = &QucsView::MouseDoNothing;
     view->MouseDoubleClickAction = &QucsView::MouseDoNothing;
     return;
     }
  }

  char *Dummy2;
  QString Dummy1;
  if(Infos) view->selComp = (*Infos) (Dummy1, Dummy2, true);

  if(view->drawn) view->viewport()->repaint();
  view->drawn = false;
  view->MouseMoveAction = &QucsView::MMoveComponent;
  view->MousePressAction = &QucsView::MPressComponent;
  view->MouseReleaseAction = &QucsView::MouseDoNothing;
  view->MouseDoubleClickAction = &QucsView::MouseDoNothing;
}

// -----------------------------------------------------------------------
// Is called when the mouse is clicked within the Content QListView.
void QucsApp::slotSelectSubcircuit(Q3ListViewItem *item)
{
  view->editText->setHidden(true); // disable text edit of component property

  if(item == 0) {   // mouse button pressed not over an item ?
    Content->clearSelection();  // deselect component in ListView
    if(view->drawn) view->viewport()->repaint();
    view->drawn = false;
    return;
  }

  if(item->parent() == 0) return;
  if(item->parent()->text(0) != tr("Schematics"))
    return;   // return, if not clicked on schematic
  if(item->text(1).isEmpty()) return;   // return, if not a subcircuit

  // delete previously selected elements
  if(view->selComp != 0)  delete view->selComp;
  if(view->selDiag != 0)  delete view->selDiag;
  view->selComp = 0;   // no component selected
  view->selDiag = 0;   // no diagram selected

  // toggle last toolbar button off
  if(activeAction) {
    activeAction->blockSignals(true); // do not call toggle slot
    activeAction->setOn(false);       // set last toolbar button off
    activeAction->blockSignals(false);
  }
  activeAction = 0;

  //view->selComp = new Subcircuit();
  view->selComp->Props.first()->Value = item->text(0);
  view->selComp->recreate();

  if(view->drawn) view->viewport()->repaint();
  view->drawn = false;
  view->MouseMoveAction = &QucsView::MMoveComponent;
  view->MousePressAction = &QucsView::MPressComponent;
  view->MouseReleaseAction = &QucsView::MouseDoNothing;
  view->MouseDoubleClickAction = &QucsView::MouseDoNothing;
}

// #######################################################################
// Is called, when "Open Project"-Button is pressed.
void QucsApp::slotProjOpenButt()
{
  view->editText->setHidden(true); // disable text edit of component property

  Q3ListBoxItem *item = Projects->selectedItem();
  if(item) slotOpenProject(item);
  else QMessageBox::information(this, tr("Info"),
				tr("No project is selected !"));
}

// #######################################################################
bool QucsApp::DeleteProject(const QString& Path, const QString& Name)
{
  view->editText->setHidden(true); // disable text edit of component property

  if(Name == ProjName) {
    QMessageBox::information(this, tr("Info"),
			tr("Cannot delete an open project !"));
    return false;
  }

  // first ask, if really delete project ?
  if(QMessageBox::warning(this, tr("Warning"),
     tr("This will destroy all the project files permanently ! Continue ?"),
     tr("&Yes"), tr("&No"), 0,1,1))  return false;

  QDir projDir = QDir(Path);
  QStringList ProjFiles = projDir.entryList("*", QDir::Files);  // all files

  // removes every file, remove("*") does not work
  QStringList::iterator it;
  for(it = ProjFiles.begin(); it != ProjFiles.end(); it++) {
     if(!projDir.remove(*it)) {
       QMessageBox::information(this, tr("Info"),
				tr("Cannot remove project file: ")+(*it));
       return false;
     }
  }

  projDir.cdUp();  // leave project directory for deleting
  if(!projDir.rmdir(Name+"_prj")) {
    QMessageBox::information(this, tr("Info"),
			     tr("Cannot remove project directory !"));
    return false;
  }

  return true;
}

// ########################################################################
// Is called, when "Delete Project"-menu is pressed.
void QucsApp::slotMenuDelProject()
{
  Q3FileDialog *d = new Q3FileDialog(QucsHomeDir.path());
  d->setCaption(tr("Choose Project Directory for Deleting"));
  d->setShowHiddenFiles(true);
  d->setMode(Q3FileDialog::DirectoryOnly);
  if(d->exec() != QDialog::Accepted) return;

  QString s = d->selectedFile();
  if(s.isEmpty()) return;

  s = s.left(s.length()-1);   // cut off trailing '/'
  int i = s.findRev('/');
  if(i > 0) s = s.mid(i+1);   // cut out the last subdirectory
  s.remove("_prj");
  DeleteProject(d->selectedFile(), s);
  readProjects();   // re-reads all projects and inserts them into the ListBox
}

// #######################################################################
// Is called, when "Delete Project"-Button is pressed.
void QucsApp::slotProjDelButt()
{
  Q3ListBoxItem *item = Projects->selectedItem();
  if(!item) {
    QMessageBox::information(this, tr("Info"),
			     tr("No project is selected !"));
    return;
  }

  if(!DeleteProject(QucsHomeDir.filePath(item->text()+"_prj"),
	item->text()))  return;
  Projects->removeItem(Projects->currentItem());  // remove from project list
}

// #######################################################################
void QucsApp::switchEditMode(bool SchematicMode)
{
  if(SchematicMode) {
    symEdit->setMenuText(tr("Edit Circuit Symbol"));
    symEdit->setStatusTip(tr("Edits the symbol for this schematic"));
    symEdit->setWhatsThis(
	tr("Edit Circuit Symbol\n\nEdits the symbol for this schematic"));
  }
  else {
    symEdit->setMenuText(tr("Edit Schematic"));
    symEdit->setStatusTip(tr("Edits the schematic"));
    symEdit->setWhatsThis(tr("Edit Schematic\n\nEdits the schematic"));
  }


  fillComboBox(SchematicMode);
  slotSetCompView(0);

  Acts.editActivate->setEnabled(SchematicMode);
//  Acts.insEquation->setEnabled(SchematicMode);
  Acts.insGround->setEnabled(SchematicMode);
//  Acts.insPort->setEnabled(SchematicMode);
  Acts.insPlot->setEnabled(SchematicMode);
  Acts.insDot_Model->setEnabled(SchematicMode);
  Acts.insDot_Top->setEnabled(SchematicMode);
  Acts.insDot_Bottom->setEnabled(SchematicMode);
  Acts.insWire->setEnabled(SchematicMode);
  Acts.insLabel->setEnabled(SchematicMode);
//  Acts.setMarker->setEnabled(SchematicMode);
  simulate->setEnabled(SchematicMode);
}

// #######################################################################
void QucsApp::changeSchematicSymbolMode(QucsDoc *d)
{
  if(!d->symbolMode) {
    switchEditMode(true);

    d->Comps  = &(d->DocComps);
    d->Wires  = &(d->DocWires);
    d->Nodes  = &(d->DocNodes);
    d->Diags  = &(d->DocDiags);
    d->Paints = &(d->DocPaints);

    QString *ps = d->UndoStack.current();
    if(ps != d->UndoStack.getFirst())  undo->setEnabled(true);
    else  undo->setEnabled(false);
    if(ps != d->UndoStack.getLast())  redo->setEnabled(true);
    else  redo->setEnabled(false);
  }
  else {
    // go into select modus to avoid placing a forbidden element
    Acts.select->setOn(true);

    switchEditMode(false);

    d->Comps  = &SymbolComps;
    d->Wires  = &SymbolWires;
    d->Nodes  = &SymbolNodes;
    d->Diags  = &SymbolDiags;
    d->Paints = &(d->SymbolPaints);

    // If the number of ports is not equal, remove or add some.
    unsigned int countPort = d->adjustPortNumbers();

    // If a symbol does not yet exist, create one.
    if(d->SymbolPaints.count() == countPort) {
      int h = 30*((countPort-1)/2) + 10;
      d->SymbolPaints.prepend(new ID_Text(-20, h+4));

      d->SymbolPaints.append(
	new GraphicLine(-20, -h, 40,  0, QPen(Qt::darkBlue,2)));
      d->SymbolPaints.append(
	new GraphicLine( 20, -h,  0,2*h, QPen(Qt::darkBlue,2)));
      d->SymbolPaints.append(
	new GraphicLine(-20,  h, 40,  0, QPen(Qt::darkBlue,2)));
      d->SymbolPaints.append(
	new GraphicLine(-20, -h,  0,2*h, QPen(Qt::darkBlue,2)));
//  Texts.append(new Text( -7,  0,"sub"));

      unsigned int i=0, y = 10-h;
      while(i<countPort) {
	i++;
	d->SymbolPaints.append(
	  new GraphicLine(-30, y, 10, 0, QPen(Qt::darkBlue,2)));
	d->SymbolPaints.at(i)->setCenter(-30,  y);

	if(i == countPort)  break;
	i++;
	d->SymbolPaints.append(
	  new GraphicLine( 20, y, 10, 0, QPen(Qt::darkBlue,2)));
	d->SymbolPaints.at(i)->setCenter(30,  y);
	y += 60;
      }
      d->sizeOfAll(d->UsedX1, d->UsedY1, d->UsedX2, d->UsedY2);
    }

    QString *ps = d->UndoSymbol.current();
    if(ps != d->UndoSymbol.getFirst())  undo->setEnabled(true);
    else  undo->setEnabled(false);
    if(ps != d->UndoSymbol.getLast())  redo->setEnabled(true);
    else  redo->setEnabled(false);
  }
}

// #######################################################################
// Is called if the "symEdit" action is activated, i.e. if the user
// switches between the two painting mode: Schematic and (subcircuit)
// symbol.
void QucsApp::slotSymbolEdit()
{
  view->editText->setHidden(true); // disable text edit of component property

  QucsDoc *d = view->Docs.current();
  d->symbolMode = !(d->symbolMode);  // change mode
  d->switchPaintMode();   // twist the view coordinates
  changeSchematicSymbolMode(d);

  // This can only be true when switching to the symbol the first time.
  if(d->UndoSymbol.isEmpty()) {
    d->setChanged(false, true); // "not changed" state, but put on undo stack
    //d->UndoSymbol.current()->at(1) = 'i';  // state of being unchanged
    //Changed by Shivam  At gives constant char in Qt4.5
    (*(d->UndoSymbol.current()))[1] = 'i';

  }

  view->resizeContents(int(d->Scale*double(d->ViewX2-d->ViewX1)),
                       int(d->Scale*double(d->ViewY2-d->ViewY1)));
  view->setContentsPos(d->PosX, d->PosY);
  view->viewport()->repaint();
  view->drawn = false;
}
