// This is an capacitor model
//
//                cap
//
//                 ||
//         o--+----||----+--o
//            |    ||    |
//            |          |
//            +--\/\/\/--+
//               int_g
//
// by Carlos E. Christoffersen

#ifndef C_h
#define C_h 1

class C : public Element
{
	public:

  C(const string& iname);

  ~C() {}

  static const char* getNetlistName()
  {
    return einfo.name;
  }

  // Do some local initialization
  virtual void init() throw(string&);

  // fill MNAM
  virtual void fillMNAM(FreqMNAM* mnam);
  virtual void fillMNAM(TimeMNAM* mnam);

  // State variable transient analysis
  virtual void svTran(TimeDomainSV *tdsv);
  virtual void deriv_svTran(TimeDomainSV *tdsv);

	private:

  // Element information
  static ItemInfo einfo;

  // Number of parameters of this element
  static const unsigned n_par;

  // Parameter variables
  double cap;
  double int_g;
  bool time_d;


  // Parameter information
  static ParmInfo pinfo[];

};

#endif

