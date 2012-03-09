// This is a Current controlled Voltage source
//
//                 2   I1              I2   0
//                  o--->----+      +--<--o
//                           |      |
//           	       	     |      |     +
//           	     	   (\|/)  (rI1)   Vout
//           	     	     |      |     -
//                           |      |
//           	    o--------+      +-----o
//                 3                        1
//
//
//the output impedance(internal resistance is assumed to be in series with the
//source and included with in.
// the polynomial functionality to enable multiple terminals
// is un-implemented

#ifndef Ccvs_h
#define Ccvs_h 1

class Ccvs : public Element
{
	public:

  Ccvs(const string& iname) ;

  ~Ccvs() {}

  static const char* getNetlistName()
  {
    return einfo.name;
  }
  
  // This element adds equations to the MNAM
  virtual unsigned getExtraRC(const unsigned& eqn_number, 
			      const MNAMType& type);
  virtual void getExtraRC(unsigned& first_eqn, unsigned& n_rows) const;
  

  // This is a multi-reference element.
  virtual void getLocalRefIdx(UnsignedVector& local_ref_vec,
	TerminalVector& term_list) ;
  // fill MNAM
  virtual void fillMNAM(FreqMNAM* mnam) ;
  virtual void fillMNAM(TimeMNAM* mnam) ;

	private:
  // Element information
  static ItemInfo einfo ;

  // Number of parameters of this element
  static const unsigned n_par ;

  // Parameter variables
  double r ;
  double ri ;
  double ro ;
  
  unsigned int my_row ;
  unsigned int n_rows ;

  // Parameter information
  static ParmInfo pinfo[] ;
} ;

#endif
