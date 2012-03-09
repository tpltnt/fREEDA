// ----------------------------------------------------------------------
// Here should be all help-related code
// ----------------------------------------------------------------------

extern "C" 
{
	#include <stdio.h>
	#include <time.h>
	#include "./inout/parser.h"
	#include "./inout/report.h"
}

#include <iostream>
#include <cstring>
#include "./network/CircuitManager.h"
#include "help.h"

using std::cout;

// Global pointer to version string
extern const char* freeda_version;

// Global pointer to the main circuit
// Used by output routines
extern Circuit* main_cir;

const char * fr_homepage = "http://www.freeda.org";

void pCatalog()
{
  ElementManager* em = ElementManager::getElementManager();
  // Print element catalog
  em->printCatalog();
}

void pCatalogElement(char *elementName)
{
  ElementManager* em = ElementManager::getElementManager();
  // Print catalog for the element elementName making sure it is lower case.
  em->printCatalogElement(elementName,1);
}


void pHelpMessage()
{
  // Print help message
  cout << "\nUsage:\n\n";
  cout << "\t\tfREEDA <file.net> [<file.out>]\n\n";
  cout << "\t\tor\n\n";
  cout << "\t\tfREEDA <flag>\n\n";
  cout << "Available Flags:\n\n";
  cout << "-a, \t\t : generate analysis html catalog file.\n";
  cout << "-c element \t : generate html catalog file of element.\n";
  cout << "-c full, --catalog full :";
  cout << " generate complete element html catalog file.\n";
  cout << "-c, --catalog : generate complete element documentaion and copy";
  cout << " pdf's from source tree to FREEDA_WEB_DOCUMENTATION/elements\n";
  cout << " caution: this can take several minutes.\n";
  cout << "-h, --help    : command line help (this message).\n";
  cout << "-l, --licence : show license information.\n";
  cout << "-s, --sanity  : do a sanity check (for developers).\n";
  cout << "-v, --version : print program version.\n\n";
  cout << "fREEDA home page: " << fr_homepage << "\n\n";
}

void pVersion()
{
  // Print version and copyright
  cout << "\n";
  cout << "-------------------------\n";
  cout << "fREEDA Circuit Simulator\n";
  cout <<"-------------------------\n";
  cout <<" Version " << freeda_version << " compiled on "
    <<__DATE__<<" "<<__TIME__<<".\n\n";
}


void pLicence()
{
  cout <<"--------------------------------------------------------------------------\n";
  cout <<"                fREEDA Circuit Simulator                                  \n";
  cout <<"--------------------------------------------------------------------------\n";
  cout <<" Copyright (C) 1998-2007 North Carolina State University \n";

  cout << " This program is free software; you can redistribute it and/or modify \n";
  cout << " it under the terms of the GNU General Public License as published by \n";
  cout << " the Free Software Foundation; either version 2 of the License, or \n";
  cout << " (at your option) any later version. \n";
  cout << "\n";
  cout << " This program is distributed in the hope that it will be useful,\n";
	cout << " but WITHOUT ANY WARRANTY; without even the implied warranty of\n";
	cout << " MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the\n";
	cout << " GNU General Public License for more details.\n";

	cout << " You should have received a copy of the GNU General Public License\n";
	cout << " along with this program; if not, write to the Free Software\n";
	cout << " Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA\n";
	cout << " fREEDA home page: http://www.freeda.org \n\n";
	cout << " fREEDA uses libraries which are distributed using their respective\n";
	cout << " licences. See each library page for more info.\n";
	cout << "\n";
	cout << " ADOL-C: http://www.math.tu-dresden.de/~adol-c\n";
	cout << " fftw: http://www.fftw.org/\n";
	cout << " NNES: http://www.netlib.org/opt/\n";
	cout << " sparse: http://www.netlib.org/sparse/\n";
	cout << "--------------------------------------------------------------------------\n";
}

void dumpOptions()
{
  // Print verbose data about options, output requests and network if
  // option verbose set.
  int  type;
  generic_t gval;

  if(St_GetSym(capOptionsT_P, "verbose", &type, &gval) == ST_SYM_FOUND) 
  {
    fprintf(output_F,
		"************************************************************************\n\n");
    fprintf(output_F, " *** Network Dump:\n");
    fprintf(output_F,"\n\n *** Title:\n\n %s\n\n",spiceTitle);
    DumpTable(output_F, capOptionsT_P);
    DumpTable(output_F, capOutputT_P);
    v_dumpExpressions();
    v_dumpSweeps();

    // Now Print circuits
    Circuit* cir;
    unsigned n = the_CM->getNumberOfCircuits();
    for (unsigned j = 0; j < n; j++) 
    {
      cir = the_CM->getCircuit(j);
      if (cir->isSubcircuit()) 
      {
				fprintf(output_F, "\n\n *** Subircuit \"%s\" listing:\n\n",
				cir->getInstanceName().c_str());
				fprintf(output_F, "     Used by %d instances\n\n",
				cir->getNumberOfUsers());
      }
      else
				fprintf(output_F, "\n\n *** Circuit \"%s\" listing:\n\n",
			cir->getInstanceName().c_str());

      Terminal* term;
      cir->setFirstElement(0);
      Element* elem = cir->nextElement();
      while(elem) 
      {
				// Print element name
				fprintf(output_F, "%s - %s\n",
				elem->getInstanceName().c_str(),
				elem->getDescription().c_str());
				// Look at the terminals of the element
				unsigned tc = elem->getNumTerms();
				for (unsigned i=0; i < tc; i++) 
        {
					term = elem->getTerminal(i);
					assert(term);
					fprintf(output_F, "            %s\n",
					term->getInstanceName().c_str());
				}
				fprintf(output_F, "\n");
				// get next element pointer
				elem = cir->nextElement();
      }
    }
    fprintf(output_F, "\n\n");
    fprintf(output_F,
		"************************************************************************\n\n");
  }
}

