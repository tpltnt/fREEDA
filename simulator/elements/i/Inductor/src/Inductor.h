// This is an inductor model
//
//                 ind    int_res
//                 
//          o-----CCCCC---\/\/\/\-----o
//
// by Carlos E. Christoffersen

#ifndef Inductor_h
#define Inductor_h 1

class Inductor : public Element
{
public:
  
  Inductor(const string& iname);

  ~Inductor() {}

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

  // State variable transient analysis
  virtual void svTran(TimeDomainSV *tdsv);
  virtual void deriv_svTran(TimeDomainSV *tdsv);

  inline double getInd()
  {
      return *((double*)(paramvalue[0]));
  }

  inline unsigned getMyrow()
  {
      return my_row;
  }

private:

  unsigned my_row;

  // Element information
  static ItemInfo einfo;

  // Number of parameters of this element
  static const unsigned n_par;
  static const double factor;

  // Parameter variables
  double ind;
  double int_res;
  bool time_d;

  // Parameter information
  static ParmInfo pinfo[];
};

#endif
