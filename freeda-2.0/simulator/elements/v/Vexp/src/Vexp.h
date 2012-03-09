// This is an Exponential voltage source
//
//       n1     + ---  -   n2
//          o----( ~ )-----o
//                ---
//
//      amplitude, frequency, phase
//
// by Satish V. Uppathil

#ifndef Vexp_h
#define Vexp_h 1

class Vexp : public Element
{
	public:

  Vexp(const string& iname) ;

  ~Vexp() {}

  static const char* getNetlistName()
  {
    return einfo.name ;
  }

  // This element adds equations to the MNAM
  virtual unsigned getExtraRC(const unsigned& eqn_number,
	const MNAMType& type) ;
  virtual void getExtraRC(unsigned& first_eqn, unsigned& n_rows) const ;

  // fill MNAM
  virtual void fillMNAM(TimeMNAM* mnam) ;
  virtual void fillSourceV(TimeMNAM* mnam) ;

  // State variable transient analysis
  virtual void svTran(TimeDomainSV *tdsv) ;
  virtual void deriv_svTran(TimeDomainSV *tdsv) ;

	private:

  // row assigned to this instance by the FreqMNAM
  unsigned my_row ;

  // Convenience variables
  double omega, phase_rad ;

  // Element information
  static ItemInfo einfo ;

  // Number of parameters of this element
  static const unsigned n_par ;

  // Parameter variables
  double v1, v2, tdr, tdf, tcr, tcf ;

  // Parameter information
  static ParmInfo pinfo[] ;

} ;

#endif

