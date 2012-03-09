#ifndef _SLIST_H_
#define _SLIST_H_

#include <iostream>
using namespace std;
#include <cstdlib>

typedef struct sllist
{
  struct sllist* next;
  double val;
  int row;
  int col;
} *Sllist;

Sllist new_sllist();
Sllist insert_sllist(int nrow,int ncol, double dval, Sllist head);
void sllist_print(Sllist head);

#endif

