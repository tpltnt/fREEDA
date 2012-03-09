// This is a pulsed voltage source with exponential taper
// (transient analysis only)
//
//       n1     + ---  -   n2
//          o----(   )-----o
//                ---
// This element behaves as a short circuit for AC/HB analysis.
//
//
// by Frank P. Hart 8/27/2003

#ifndef Vpulsexp_h
#define Vpulsexp_h 1

class Vpulsexp : public Element
{
	public:

  Vpulsexp(const string& iname);

  ~Vpulsexp() {}

  static const char* getNetlistName()
  {
    return einfo.name;
  }

  // This element adds equations to the MNAM
  virtual unsigned getExtraRC(const unsigned& eqn_number,
	const MNAMType& type);
  virtual void getExtraRC(unsigned& first_eqn, unsigned& n_rows) const;

  // fill MNAM
  virtual void fillMNAM(FreqMNAM* mnam);
  virtual void fillMNAM(TimeMNAM* mnam);
  virtual void fillSourceV(TimeMNAM* mnam);

  // State variable transient analysis
  virtual void svTran(TimeDomainSV *tdsv);
  virtual void deriv_svTran(TimeDomainSV *tdsv);

	private:

  // row assigned to this instance by the FreqMNAM
  unsigned my_row;

  // Time for the integer number of periods before current time
  double int_per_time;

  // Element information
  static ItemInfo einfo;

  // Number of parameters of this element
  static const unsigned n_par;

  // Parameter variables
  double v1, v2, td, tr, tf, pw, per;

  // Parameter information
  static ParmInfo pinfo[];

  // Time constants derived from tr and tf
  // Normalized non-zero rail voltage vn.
  double tau_r, tau_f, vn;

  // Pulse width high and low vars for better readability
  double pwl, pwh;

};

#endif

