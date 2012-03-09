// This may look like C code, but it is really -*- C++ -*-
//
// This is a Amplitude Modulated current source 
//
//       n1       ---      n2
//          o----( < )-----o
//                ---     
//
// by Satish V. Uppathil

#ifndef IAM_h
#define IAM_h 1

class IAM : public Element
{
public:
  
  IAM(const string& iname);

  ~IAM() {}

  static const char* getNetlistName();

  // Do some local initialization
  //virtual void init() throw(string&);

  // fill MNAM
  //virtual void fillMNAM(FreqMNAM* mnam);
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
  double oc, sa, fcarrier, fmod, td ;
  
  // Parameter information
  static ParmInfo pinfo[];

};
#endif


