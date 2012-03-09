// This may look like C code, but it is really -*- C++ -*-
//
// This is a Single Frequency FM current source 
//
//       n1       ---      n2
//          o----( < )-----o
//                ---     
//
//      amplitude, frequency, phase
//
// by Satish V. Uppathil

#ifndef IsfFM_h
#define IsfFM_h 1

class IsfFM : public Element
{
public:
  
  IsfFM(const string& iname);

  ~IsfFM() {}

  static const char* getNetlistName()
  {
    return einfo.name;
  }

  // Do some local initialization
  // virtual void init() throw(string&);

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

  // Number of parameters of this element
  static const unsigned n_par;

  // Parameter variables
  double io, ia, fcarrier, mdi, fsignal ;

 // Parameter information
  static ParmInfo pinfo[];

};
#endif


