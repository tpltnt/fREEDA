// This is a DC and sinusoidal voltage source 
//
//       n1     + ---  -   n2
//          o----( ~ )-----o
//                ---     
//
//      amplitude, frequency, phase
//
// by Carlos E. Christoffersen

#ifndef Vsource_h
#define Vsource_h 1

class Vsource : public Element
{
public:
  
  Vsource(const string& iname);

  ~Vsource() {}

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
  double omega, phase_rad, ramp_slope;

  // Element information
  static ItemInfo einfo;

  // Number of parameters of this element
  static const unsigned n_par;

  // Parameter variables
  double vdc, vac, frequency, phase, delay, tr;
  double periods;

  // Parameter information
  static ParmInfo pinfo[];

};

#endif
