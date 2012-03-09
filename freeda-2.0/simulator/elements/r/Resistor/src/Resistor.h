// This is an resistor model
//
//                res
//         o----/\/\/\----o
//
// by Carlos E. Christoffersen

#ifndef Resistor_h
#define Resistor_h 1

class Resistor : public Element
{
	public:
	
  Resistor(const string& iname);
  ~Resistor() {}
	
  static const char* getNetlistName()
  {
    return einfo.name;
  }
	
  // fill MNAM
  virtual void fillMNAM(FreqMNAM* mnam);
  virtual void fillMNAM(TimeMNAM* mnam);
	
	private:
	
  // Element information
  static ItemInfo einfo;
	
  // Number of parameters of this element
  static const unsigned n_par;
	
  // Parameter variables
  double res;
	
  // Parameter information
  static ParmInfo pinfo[];
};

#endif

