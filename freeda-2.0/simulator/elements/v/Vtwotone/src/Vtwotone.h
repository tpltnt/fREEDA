// This is a twotone sinusoidal voltage source 
//
//       n1     + ---  -   n2
//          o----( ~ )-----o
//                ---     
//
//      amplitude, frequency, phase
//
// by Aravind K. Mikkilineni

#ifndef Vtwotone_h
#define Vtwotone_h 1

class Vtwotone : public Element
{
public:
  
  Vtwotone(const string& iname);

  ~Vtwotone() {}

  static const char* getNetlistName()
  {
    return einfo.name;
  }

  // Do some local initialization
  virtual void init() throw(string&);

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

  // Convenience variables
  double omega1, omega2, phase_rad1, phase_rad2;

  // Element information
  static ItemInfo einfo;

  // Number of parameters of this element
  static const unsigned n_par;

  // Parameter variables
  double vac1, f1, phase1, vac2, f2, phase2, delay;

  // Parameter information
  static ParmInfo pinfo[];

};

#endif
