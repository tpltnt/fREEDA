// Class definition for the AC analysis
// Author:
// Carlos E. Christoffersen
// logsweep functionality contributed by Kris Andersson

#ifndef AC_h
#define AC_h 1

#include "Analysis.h"

// Main class definition follows

class AC : public Analysis
{
public:

  AC();

  ~AC();

  // The main analysis routine.
  virtual void run(Circuit* cir);
  
  //Return the name of the element
  static const char* getNetlistName()
  {
     return ainfo.name;
  }

private:

  // Write results to network.
  void doOutput();

  // ------------------------------  Variables

  // Circuit pointer
  Circuit* my_cir;
  // Frequency vector
  DenseDoubleVector f_vec;
  // Result vector of vectors
  DenseComplexVector* result_vec;

  // ------------------- Parameter-related variables

  // Analysis information
  static ItemInfo ainfo;

  // Number of parameters of this element
  static const unsigned n_par;

  // Analysis parameters
  double f_start;
  double f_stop;
  int n_freqs;
  bool logsweep; // thanks to kris andersson

  // Parameter information
  static ParmInfo pinfo[];
};

#endif

