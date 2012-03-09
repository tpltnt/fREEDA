//---------------------------------------------------------------------
// This is an analog nth-order type-2 chebyshev bandstop filter model
//
//         --------------
//         |   n,   |
//  n1 o---|   fc, bw,  |---o n3
//         | z0, ripple |
//         |            |
//         --------------
//                |
//                |
//                o n2
//
// Michael Steer, September 2008
//
// This model designs a type-2 Cauer topology chebyshev bandstop
// filter, using parallel inductors and capacitors and series
// inductors and capacitors in a ladder structure.  The filter can be
// of any order up to a maximum order of 100.  The filter has a center
// frequency, f0, a bandwidth, bw, and a passband ripple, ripple. The
// filter is designed so that the input and output impedance, z0, of
// the filter is the same.  
//
// The native fREEDA netlist format is
// Chebyshevbsf:1 n1 n2 n3 n=? f0=? bw=? z0=? ripple=? q=? 
// where:
// Chebyshevbsf indicates the model to use
// Chebyshevbsf:1 is the name of this instance
// n1 is the input terminal (Can be a string or an integer)
// n2 is the reference terminal (Can be a string or an integer)
// n3 is the output terminal (Can be a string or an integer)
// n is the filter order
// f0 is the center frequency of the filter in Hz
// bw is the bandwidth of the filter in Hz
// z0 is the input/output impedance in Ohms
// ripple is the passband ripple in decibels
// q is the Q of the resonators
//----------------------------------------------------------------------

// This file sets up the class for the Chebyshevbsf element

// The following is done in case Chebyshevbsf.h is included twice. If
#ifndef ChebyshevBSF_h
#define ChebyshevBSF_h 1

class ChebyshevBSF : public Element
{
	public:
// This is the prototype of the creator routine
  ChebyshevBSF(const string& iname);

// This is the prototype of the destructor routine
  ~ChebyshevBSF() {}

  static const char* getNetlistName()
  {
    return einfo.name;
  }

  // local initialization to calculate L's and C's and to build up the
  // filter
  virtual void init() throw(string&);

  // Get a vector with the indexes of the local reference nodes.
  virtual void getLocalRefIdx(UnsignedVector& local_ref_vec,
			      TerminalVector& term_list);

  //  For the frequency domain form (complex elements)
  virtual void fillMNAM(FreqMNAM* mnam);

  //  For the time domain form (real elements)
  virtual void fillMNAM(TimeMNAM* mnam);

	private:

  // Set up the variable to store element information
  static ItemInfo einfo;

  // Set up the variable to store the number of parameters of this element
  static const unsigned n_par;

  // List of parameters
  int n;
  double f0, bw, z0, q, ripple; 
  
  // Parameter information. Space is allocated for the pointer to the
  // pinfo vector.
  static ParmInfo pinfo[];
};

#endif

