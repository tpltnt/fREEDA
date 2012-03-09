//---------------------------------------------------------------------
// This is an analog nth-order type-2 butterworth bandpass filter model
//
//         --------------
//         |   order,   |
//  n1 o---|   fc, bw,  |---o n2
//         |     z0     |
//         |            |
//         --------------
//                |
//                |
//                o n3
//
// by Shawn D. Evans
//
// This model designs a type-2 Cauer topology butterworth bandpass
// filter, using parallel inductors and capacitors and series
// inductors and capacitors in a ladder structure.  The filter can be
// of any order up to a maximum order of 100.  The filter has a center
// frequency, fc, and a bandwidth, bw. The filter is designed so that
// the input and output impedance, z0, of the filter is the same.  
//
// The native fREEDA netlist format is
// Butterworthbpf:1 n1 n2 n3 order=? fc=? bw=? z0=?
// where:
// Butterworthbpf indicates the model to use
// Butterworthbpf:1 is the name of this instance
// n1 is the input terminal (Can be a string or an integer)
// n2 is the output terminal (Can be a string or an integer)
// n3 is the reference terminal (Can be a string or an integer)
// order is the filter order
// fc is the center frequency of the filter in Hz
// bw is the bandwidth of the filter in Hz
// z0 is the input/output impedance in Ohms
//----------------------------------------------------------------------

// This file sets up the class for the Butterworthbpf element

// The following is done in case ButterworthBPF.h is included twice. If
#ifndef ButterworthBPF_h
#define ButterworthBPF_h 1

class ButterworthBPF : public Element
{
	public:
// This is the prototype of the creator routine
  ButterworthBPF(const string& iname);

// This is the prototype of the destructor routine
  ~ButterworthBPF() {}

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
  int order;
  double fc, bw, z0; 
  
  // List of private variables.
  double wc, bw_rad, g[100], C[100], L[100];

  // Parameter information. Space is allocated for the pointer to the
  // pinfo vector.
  static ParmInfo pinfo[];
};

#endif

