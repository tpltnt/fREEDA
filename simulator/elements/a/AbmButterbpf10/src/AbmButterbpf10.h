// This is for a 4-terminal device defined by y = H(s) = I(s) / V(s)
//
//           +------+
//  1 o------+      +------o 2
//     out   +      +control
//  0 o------+      +------o 3
//	     +------+
//

#ifndef AbmButterbpf10_h
#define AbmButterbpf10_h 1

class AbmButterbpf10 : public Element
{
	public:

  AbmButterbpf10(const string& iname);

  ~AbmButterbpf10() {};

  static const char* getNetlistName()
  {
    return einfo.name;
  }

  //Do some local initialization
  virtual void init() throw(string&);

  virtual unsigned getExtraRC(const unsigned& eqn_number, const MNAMType& type);

  virtual void getExtraRC(unsigned& first_eqn, unsigned& n_rows) const;

  // fill MNAM
  virtual void fillMNAM(FreqMNAM* mnam);
  virtual void fillMNAM(TimeMNAM* mnam);

	private:

  // Element information
  static ItemInfo einfo;

  // Number of parameters of this element
  static const unsigned n_par;

  // Parameter variables
  double bw, f0;

  // Internal variables
  unsigned my_start_row, extra_rcs;
  double bw_rad, w0, w02, w04;
  DenseDoubleVector test;

  DenseDoubleVector a, b, c;	// Final denom polynomial coefficients
  DenseDoubleVector d, e, f, h;
  DenseDoubleVector k0;
  DenseDoubleVector k1;			// Prototype num polynomial coefficients
  DenseDoubleVector bt;			// Prototype den polynomial coefficients

  // Parameter information
  static ParmInfo pinfo[];
};

#endif

