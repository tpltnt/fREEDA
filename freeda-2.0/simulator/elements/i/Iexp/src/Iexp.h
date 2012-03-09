// This may look like C code, but it is really -*- C++ -*-
//
// This is a Exponential current source 
//
//       n1       ---      n2
//          o----( < )-----o
//                ---     
//      amplitude, frequency, phase
//
// by Satish V. Uppathil

#ifndef Iexp_h
#define Iexp_h 1

class Iexp : public Element
{
public:
  
  Iexp(const string& iname);

  ~Iexp() {}

  static const char* getNetlistName()
  {
    return einfo.name;
  }

  // Do some local initialization
  //virtual void init() throw(string&);

  // fill MNAM
  virtual void fillMNAM(FreqMNAM* mnam);
  virtual void fillMNAM(TimeMNAM* mnam);
  virtual void fillSourceV(TimeMNAM* mnam);

  // State variable transient analysis
  virtual void svTran(TimeDomainSV *tdsv);
  virtual void deriv_svTran(TimeDomainSV *tdsv);

private:

  // Element information
  static ItemInfo einfo;

  // Convenience variables
  double omega, phase_rad ;

  // Number of parameters of this element
  static const unsigned n_par;

  // Parameter variables
  double v1, v2, tdr, tdf, tcr, tcf ;


 // Parameter information
  static ParmInfo pinfo[];

};
#endif


