// This is for a N-terminal device defined by y = H(s) = I(s) / V(s)
//
//  		     +----------+
//      2 o----------+          +---------o 3
//  		     +          +
//      1 o----------+          +---------o 4
//  		     +          +
//      0 o----------+          +---------o 5
//  	  	|    +          +    |
//  	  	|    +          +    |
//              |    +          +    |
//              |    +          +    |
//              |    +          +    |
//              |    +          +    |
//  		|    +----------+    |
//  		          |
//  		          |
//  		          o
//  		          Nth reference terminal

#ifndef FosterNPort_h
#define FosterNPort_h 1

class FosterNPort : public Element
{
  public:

  FosterNPort(const string& iname);

  ~FosterNPort() {};

  static const char* getNetlistName()
  {
    return einfo.name;
  }

  //Do some local initialization

  double Value(char * line);

  virtual void init() throw(string&);

  virtual unsigned getExtraRC(const unsigned& eqn_number, const MNAMType& type);

  virtual void getExtraRC(unsigned& first_eqn, unsigned& n_rows) const;

  // fill MNAM
  virtual void fillMNAM(FreqMNAM* mnam);
  virtual void fillMNAM(TimeMNAM* mnam);

  private:

  //Work variables

  // Number of complex poles and real poles for each element
  int *noCpoles; // number if complex poles for element
  double **cPoleRe; //  Pointer to Pointer of real part of complex poles
  double **cPoleIm; //  Pointer to Pointer of imaginary part of complex poles
  double **cResidueRe; //  Pointer to Pointer of real part of complex residues
  double **cResidueIm; //  Pointer to Pointer of imaginary part of residues
  int *noRpoles; // number of real polese for element
  double **rPole; //  Pointer to Pointer of real part of complex poles
  double **rResidue; //  Pointer to Pointer of real part of complex residues

  // Element information
  static ItemInfo einfo;

  // Number of parameters of this element
  static const unsigned n_par;

  // Parameter variables
  string filename;
  int ports;

  // indeces tracking MNAM entries
  unsigned my_start_row, extra_rcs;
  // unsigned extra_rcs, i_start_row;

  // Parameter information
  static ParmInfo pinfo[];

};

#endif

