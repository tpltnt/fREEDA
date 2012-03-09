// This is a Current-controlled Current source
//
//                 2                        0
//                  o--------+      +------------o
//                           |      |          
//           	      	       |      |           
//           	       Iin    \|/    \|/ Iout    ro
//           	    	         |      |          
//                           |      |          
//           	      o--------+      +------------o
//                 3                        1
//
//
//the output impedance of the Current-controlled current source
//is assumed to be in parallel with the output current source

#ifndef Cccs_h
#define Cccs_h 1

class Cccs : public Element
{
	public:

  Cccs(const string& iname) ;

  ~Cccs() {}

  static const char* getNetlistName()
  {
    return einfo.name;
  }

  // This is a multi-reference element.
  virtual void getLocalRefIdx(UnsignedVector& local_ref_vec,
	TerminalVector& term_list) ;
	
  virtual unsigned getExtraRC(const unsigned& eqn_number, 
                              const MNAMType& type);
  virtual void getExtraRC(unsigned& first_eqn, unsigned& n_rows) const;	
	
	
  // fill MNAM
  virtual void fillMNAM(FreqMNAM* mnam) ;
  virtual void fillMNAM(TimeMNAM* mnam) ;

	private:
  // Element information
  static ItemInfo einfo ;

  // Number of parameters of this element
  static const unsigned n_par ;

  // Parameter variables
  double a ;
  double ro ;
  
  unsigned int my_row ;
  unsigned int n_rows ;

  // Parameter information
  static ParmInfo pinfo[] ;
} ;

#endif
