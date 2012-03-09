//Interface class


#ifndef INTERFACE_H
#define INTERFACE_H

#include <qmessagebox.h>
#include <qstring.h>
#include <string>
using namespace std;

#include "components/component.h"

extern int maimed;

class Interface : public Component  {
public: 
  Interface(QString& elem_type);
  ~Interface();

  Component* newOne();
  static Component* info(QString&, char* &, bool getNewOne=false);
  void drawGenericElement();
  void drawGenericAnalysis();

private:
  
  int numberOfTerminals;
  unsigned numOfPrms; //Number of Parameters  
  bool display;
  string  parname, comment, type, dflt_val, required;
  //string nport_str; //Work around for nports having dynamic ports
  // changed by shivam  Changed elem_type & nport_str from string to QString
  QString elem_type,nport_str;

};

#endif
