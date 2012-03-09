// This is a Voltage controlled current source
//
//                 3                        0
//                  o--------+      +---+--o
//                           |      |   |
//           	    +  	     \      |   \      +
//           	  Vin	     / Ri  ( )  / Ro   Vout
//           	    -	     \      |   \      -
//                           |      |   |
//           	    o--------+      +---+--o
//                  2                       1
//
// by Carlos E. Christoffersen
// specify the output terminals before the input terminals
// in the netlist (NMK)

#ifndef Vccs_h
#define Vccs_h 1

class Vccs : public Element
{
	public:

  Vccs(const string& iname);

  ~Vccs() {}

  static const char* getNetlistName()
  {
    return einfo.name;
  }

  // This is a multi-reference element.
  virtual void getLocalRefIdx(UnsignedVector& local_ref_vec,
	TerminalVector& term_list);
  // fill MNAM
  virtual void fillMNAM(FreqMNAM* mnam);
  virtual void fillMNAM(TimeMNAM* mnam);

	private:

  // Element information
  static ItemInfo einfo;

  // Number of parameters of this element
  static const unsigned n_par;

  // Parameter variables
  double g;
  double ri;
  double ro;

  // Parameter information
  static ParmInfo pinfo[];

};

#endif

