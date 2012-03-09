#include "AC.h"
#include "FreqMNAM.h"

extern "C"
{
  #include "../inout/ftvec.h"
  #include "../inout/report.h"
}

#include <cstdio>

// Static members
const unsigned AC::n_par = 4;

// Element information
ItemInfo AC::ainfo =
{
  "ac",
  "AC Analysis",
  "Carlos E. Christoffersen, Kristoffer Andersson",
  DEFAULT_ADDRESS
};

// Parameter information
ParmInfo AC::pinfo[] =
{
  {"start", "Start frequency (Hz)", TR_DOUBLE, true},
  {"stop", "Stop frequency (Hz)", TR_DOUBLE, true},
  {"n_freqs", "Number of frequencies", TR_INT, true},
  {"logsweep", "Logarithmic sweep", TR_BOOLEAN, false}
};


AC::AC() : Analysis(&ainfo,pinfo, n_par)
{
  // Parameter stuff
  paramvalue[0] = &f_start;
  paramvalue[1] = &f_stop;
  paramvalue[2] = &n_freqs;
  paramvalue[3] = &logsweep;
  result_vec = NULL;
}

AC::~AC()
{
  // For now it seems that nothing needs to be done here.
}

void AC::run(Circuit* cir)
{
  my_cir = cir;
  // Char array to store output messages
  char msg[80];
  report(MESSAGE, "*** AC Analysis ***");
  sepLine();
  newLine();

  // ----------------------------------------------------------------------
  // -------------- Set up and check parameters
  // ----------------------------------------------------------------------

  checkParams();

  // ----------------------------------------------------------------------
  // -------------- Set up frequency vector.
  // ----------------------------------------------------------------------
  // Create frequency vector
  // First element is f_start
  // Last element is f_stop
  if(logsweep)
  {
    f_vec.resize(n_freqs);
    for (int i = 0; i < (n_freqs - 1); i++)
      f_vec[i] = pow(10.0,f_start + i * (f_stop - f_start)/(n_freqs - 1.0));
    f_vec[n_freqs - 1] = pow(10.0, f_stop);
    sprintf(msg, "Logarithmic frequency sweep");
  }
  else
  {
    double fstep = (f_stop - f_start) / (n_freqs - 1);
    f_vec.resize(n_freqs);
    for (int i = 0; i < n_freqs; i++)
      f_vec[i] = f_start + fstep * i;

    sprintf(msg, "Frequency step = %g Hz", fstep);
  }
  report(MESSAGE, msg);

  // ----------------------------------------------------------------------
  // -------------- Build MNAM and solve at each frequency
  // ----------------------------------------------------------------------

  // Set flags to linear elements for now. Later add treatment for
  // small-signal models.
  ElemFlag mask = LINEAR;

  // Right now just build all the matrices at the same time.
  // Later we may want to build just one at a time as in
  // transient analysis.
  FreqMNAM* mnam = new FreqMNAM(f_vec, my_cir, mask);
  
  // this step is required
  mnam->FillComplete();

  // Check if there is no memory allocated for result_vec from a
  // previous run.
  if (result_vec)
    delete [] result_vec;
	result_vec = new DenseComplexVector[n_freqs];
	for (int i = 0; i < n_freqs; i++)
		result_vec[i] = DenseComplexVector(mnam->getDim());

  // Loop through all the frequencies
  for (unsigned findex=0; findex < unsigned(n_freqs); findex++)
	{
    // Solve for default source vector.
    mnam->solve(findex, result_vec[findex]);
  }

  // Delete mnam since it is no longer necessary
  delete mnam;

  // Write output
  doOutput();

  // Erase the result vector
  delete [] result_vec;
}

void AC::doOutput()
{
  // First check if the result matrices contain any data
  assert(result_vec);

  report(MESSAGE, "--- Writing output vectors ...");

  // Allocate and copy global freq vector
  allocFreqV_P(n_freqs);
  for (int findex = 0; findex < n_freqs; findex++)
    FreqV_P[findex] = f_vec[findex];

  // Temporary vector for voltages (initialized to zero)
  DenseComplexVector tmp_v(n_freqs);

  // For each terminal, assign voltage vector
  Terminal* term = NULL;
  my_cir->setFirstTerminal();
  while((term = my_cir->nextTerminal()))
	{
    // Get MNAM index
    if (term->getRC())
		{
      // Remember that MNAM indices begin at 0
      int i = term->getRC() - 1;
      // Get vector from VI
      for (int findex=0; findex < n_freqs; findex++)
        tmp_v[findex] = result_vec[findex][i];
    }
    else // This is a reference terminal
      tmp_v.putScalar(double_complex(zero));

    // Set terminal vector
    term->getTermData()->setPhasorV(tmp_v);
  }

  // Temporary vectors (initialized to zero)
  DenseComplexVector tmp_i(n_freqs);

  ElemFlag mask = LINEAR;
  // Loop throw all selected elements in circuit and fill current
  // vector if needed.
  my_cir->setFirstElement(mask);
  Element* elem = my_cir->nextElement();
  unsigned first_eqn, no_eqn;
  while(elem)
	{
    elem->getExtraRC(first_eqn, no_eqn);
    if (first_eqn)
		{
      // Fill current vector(s) in element
      for (unsigned i = 0; i < no_eqn; i++)
			{
	       // Remember that MNAM indices begin at 0
         int row = first_eqn - 1 + i;
         for (int findex = 0; findex < n_freqs; findex++)
           tmp_i[findex] = result_vec[findex][row];
         elem->getElemData()->setPhasorI(i, tmp_i);
      }
    }
    // get next element pointer
    elem = my_cir->nextElement();
  }
}

