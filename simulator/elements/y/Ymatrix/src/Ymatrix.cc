#include "../../../../network/ElementManager.h"
#include "../../../../analysis/FreqMNAM.h"
#include "Ymatrix.h"

#include <cstdio>
#include <cstdlib>
extern "C" {
#include "../../../../compat/in_interp.h"
	   }
#include "../../../lex.Ymatrix_.c"

// Static members
const unsigned Ymatrix::n_par = 2;

// Element information
ItemInfo Ymatrix::einfo = {
  "ymatrix",
  "Y-matrix defined element - Obsolete: Use nport instead",
  "Carlos E. Christoffersen, Mark Summers",
  DEFAULT_ADDRESS"category:multiport",
  "2000_08_20"
};

// Parameter information
ParmInfo Ymatrix::pinfo[] = {
  {"filename",
   "File containing the nodal-based y parameter matrix.", TR_STRING, true},
  {"pfilename",
   "Name of the file to generate with the port-based y parameter matrix.",
   TR_STRING, false}
};


Ymatrix::Ymatrix(const string& iname) : Element(&einfo, pinfo, n_par, iname)
{
  // Standard construction code.
  paramvalue[0] = &filename;
  paramvalue[1] = &(pfilename);


  setFlags(LINEAR | MULTI_REF | TR_FREQ_DOMAIN);

  // Set initial dimension for internal frequency vector.
  vec_capacity = 300;
  vec_size = 0;

  freq_vec = NULL;
  ymat_vec = NULL;

  // Set initial terminal and reference count from file.
  terms = 0;
  current_group = 0;
  // Set last_index to zero initially
  last_index = 0;
}


Ymatrix::~Ymatrix()
{
  if (freq_vec) {
    // Free allocated matrices
    for (unsigned i=0; i < vec_size; i++)
      delete ymat_vec[i];

    free(ymat_vec);
    free(freq_vec);
  }
}


void Ymatrix::check() throw(string&)
{
  // We need to check parameters first, before attepmting to open
  // the data file.
  try {
    checkParams();
  }
  catch(string& msg) {
    throw(getInstanceName() + ": " + msg);
  }
}


void Ymatrix::init() throw(string&)
{
  // Now allocate memory for the vectors. Use malloc() to later use
  // realloc()
  freq_vec = (double*)(malloc(sizeof(double) * vec_capacity));
  ymat_vec = (ComplexDenseMatrix**)(malloc(sizeof(ComplexDenseMatrix*) * vec_capacity));

  FILE* fp;

  // Open input file
  if( (fp = fopen(filename.c_str(), "r")) == NULL)
    throw(getInstanceName() + ": Cannot open input file \"" +
	  filename +"\".");

  // Point scanner pointer here.
  scan_YM.thisYM = this;
  Ymatrix_in = fp;
  // Set the line counter to one.
  scan_YM.linecount = 1;

  // Parse first part (ports) of input file
  if(Ymatrix_lex())
    throw(getInstanceName() + ": in file \"" + filename + "\": " +
	  + scan_YM.serror);

  // Check that the number of terminals in netlist is consistent
  // with number of terminals in file
  if (getTermCount() != terms + group_vec.size())
    throw(string(getInstanceName() + ": Incorrect number of terminals."));
  // If control reaches here, then the number of terminals in
  // netlist is correct.
  setNumTerms(getTermCount());

  // Continue parsing the rest of the file (y parameters).
  if (Ymatrix_lex())
    throw(getInstanceName() + ": in file \"" + filename + "\": " +
	  + scan_YM.serror);

  // Close input file
  fclose(fp);

  if (isSet(&pfilename)) {

    // Add some code here to generate a port-based version of the file
    fp = fopen(pfilename.c_str(), "w");
    for (unsigned findex = 0; findex < vec_size; findex++) {
      fprintf(fp, "%g\n", freq_vec[findex]);
      for (unsigned row = 0; row < terms; row++)
	for (unsigned col = 0; col < terms; col++)
	  fprintf(fp, "%g  %g\n", (*ymat_vec[findex])(row, col).real(),
		  (*ymat_vec[findex])(row, col).imag());
      fprintf(fp, "\n");
    }
    fclose(fp);
  }

//   // Print the ports and references to stdout
//   cout << "\n " << getInstanceName() << " terminal information: \n\n";
//   UnsignedVector lrv;
//   TerminalVector tv;
//   getLocalRefIdx(lrv, tv);
//   unsigned ngroups = lrv.size();
//   unsigned j = 0;
//   for (unsigned i=0; i < ngroups; i++) {
//     while (lrv[i] != j) {
//       cout << "    \"" << tv[j]->getInstanceName() << "\"\n";
//       j++;
//     }
//     cout << "    Local reference: \"" <<
//       tv[j]->getInstanceName() << "\" \n" << endl;
//     j++;
//   }

}

void Ymatrix::addPort(const unsigned& group_number)
{
  if (!terms || current_group != group_number) {
    // This is the first terminal or
    // there is a new group.
    // Add a new group with one count (the current terminal).
    group_vec.push_back(1);
    current_group = group_number;
  }
  else
    // Increment the terminal counter for this group.
    (group_vec.back())++;

  // Increment the number of terminals given in the file.
  terms++;
}


ComplexDenseMatrix* Ymatrix::newMatrix(const double& frequency)
{
  if (vec_size == vec_capacity) {
    // Allocate more space
    vec_capacity *= 2;
    freq_vec = (double*)(realloc(freq_vec,
				 sizeof(double) * vec_capacity));
    ymat_vec = (ComplexDenseMatrix**)
      (realloc(ymat_vec, sizeof(ComplexDenseMatrix*) * vec_capacity));
  }
  freq_vec[vec_size] = frequency;
  ymat_vec[vec_size] = new ComplexDenseMatrix(getNumTerms(), getNumTerms());
  vec_size++;
  return ymat_vec[vec_size - 1];
}


void Ymatrix::getLocalRefIdx(UnsignedVector& local_ref_vec,
			     TerminalVector& term_list)
{
  // Make sure the vectors are empty
  term_list.erase(term_list.begin(), term_list.end());
  local_ref_vec.erase(local_ref_vec.begin(), local_ref_vec.end());

  unsigned termcount = 0;
  // Iterate for each group
  for (unsigned i = 0; i < group_vec.size(); i++) {
    // Iterate for each terminal in group
    for (unsigned j = 0; j < group_vec[i]; j++) {
      term_list.push_back(getTerminal(termcount));
      // Go to next terminal (in netlist connection).
      termcount++;
    }
    // Add now the local reference terminal
    term_list.push_back(getTerminal(terms + i));
    local_ref_vec.push_back(termcount + i); // Local reference index
  }
}


void Ymatrix::fillMNAM(FreqMNAM* mnam)
{
  const double& freq(mnam->getFreq());
  unsigned i,j;
  // Find the frequency index.
  if (vec_size == 1)
    last_index = 0;
  else
    doubleHunt(freq_vec, int(vec_size), freq, &last_index);

  if (unsigned(last_index) >= vec_size) {
    if (last_index)
      // Desired frequency is too high
      // Use last information available.
      // We should also generate a warning message.
      last_index = vec_size - 1;

    for (i=0; i<getNumTerms(); i++)
      for (j=0; j<getNumTerms(); j++)
	mnam->setElement(getTerminal(i)->getRC(),
			 getTerminal(j)->getRC(),
			 (*ymat_vec[last_index])(i,j));
  }
  else if (last_index < 1) {
    // Use initial data
    for (i=0; i<getNumTerms(); i++)
      for (j=0; j<getNumTerms(); j++)
	mnam->setElement(getTerminal(i)->getRC(),
			 getTerminal(j)->getRC(),
			 (*ymat_vec[0])(i,j));
  }
  else if (freq == freq_vec[last_index])
    // Desired frequency equals one of the data frequencies
    // Set element in MNAM
    for (i=0; i<getNumTerms(); i++)
      for (j=0; j<getNumTerms(); j++)
	mnam->setElement(getTerminal(i)->getRC(),
			 getTerminal(j)->getRC(),
			 (*ymat_vec[last_index])(i,j));
  else
    // Interpolate between last_index and last_index-1
    for (i=0; i<getNumTerms(); i++)
      for (j=0; j<getNumTerms(); j++) {

	// Calculate y parameter using interpolation
	double_complex yval((*ymat_vec[last_index])(i,j));
	yval -= (*ymat_vec[last_index-1])(i,j);
	yval *= ((freq - freq_vec[last_index-1]) /
	  (freq_vec[last_index] - freq_vec[last_index-1]));
	yval += (*ymat_vec[last_index-1])(i,j);

	// Set element
	mnam->setElement(getTerminal(i)->getRC(),
			 getTerminal(j)->getRC(),
			 yval);
      }

}


