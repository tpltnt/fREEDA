// This is an resistor model
//
//                R
//         o----/\/\/\----o
//
// by Carlos E. Christoffersen
// This file sets up the class for the R element

// The following is done in case R.h is included twice.
#ifndef R_h
#define R_h 1

class R : public Element
{
	public:

  // This is the prototype of the constructor routine
  R(const string& iname);

  // This is the prototype of the destructor routine
  ~R() {}

  static const char* getNetlistName()
  {
    return einfo.name;
  }

  // Function overloading is used to fil the MNAM
  // (Modified Nodal Admittance Matrix). There are two routines here
  // of identical names. Which is called depends on the argument type
  // fill the frequency domain form of MNAM

  // For the frequency domain form (complex elements)
  virtual void fillMNAM(FreqMNAM* mnam);

  // For the tiem domain form (real elements)
  virtual void fillMNAM(TimeMNAM* mnam);

	private:

  // Set up the variable to store element information
  static ItemInfo einfo;

  // Set up the variable to store the number of parameters of this element
  static const unsigned n_par;

  // List of parameter (there is only one for this element).
  // Declare each of the parameters of an element here.
  double res;

  // Parameter information
  // Space is allocated for the pointer to the pinfo vector
  static ParmInfo pinfo[];
};

#endif

