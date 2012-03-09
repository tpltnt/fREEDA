// This is a current source
//
//       n1       ---      n2
//          o----( < )-----o
//                ---
//
//      amplitude, frequency, phase
//
// by Carlos E. Christoffersen

#ifndef Isource_h
#define Isource_h 1

class Isource : public Element
{
	public:

  Isource(const string& iname);

  ~Isource() {}

  static const char* getNetlistName()
  {
    return einfo.name;
  }

  // Do some local initialization
  virtual void init() throw(string&);

  // fill MNAM
  virtual void fillMNAM(FreqMNAM* mnam);
  virtual void fillMNAM(TimeMNAM* mnam);
  virtual void fillSourceV(TimeMNAM* mnam);

  // State variable transient analysis
  virtual void svTran(TimeDomainSV *tdsv);
  virtual void deriv_svTran(TimeDomainSV *tdsv);

	private:

  // Convenience variables
  double omega, phase_rad, ramp_slope;

  // Element information
  static ItemInfo einfo;

  // Number of parameters of this element
  static const unsigned n_par;

  // Parameter variables
  double idc, iac, frequency, phase, delay, tr, periods;

  // Parameter information
  static ParmInfo pinfo[];

};

#endif

