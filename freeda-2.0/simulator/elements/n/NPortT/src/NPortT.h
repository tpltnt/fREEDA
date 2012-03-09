// NPortT: linear multiport element given by its y parameter transfer
// function matrix.  The number of groups and terminals in each group
// is variable.
//
// Author: Carlos E. Christoffersen
//
//
//              +-----------+
//         o----|           |
// group1  o----|           |----o  group3
//         o----|           |----o
//         o----|  NPortT    |----o
//         |    |           |    |
//         ^    |           |    ^ local reference
//              |           |
// group2  o----|           |----o  ...
//         o----|           |----o  group(m)
//         |    +-----------+    |
//         ^                     ^
//
// Terminal order is: Depends on the order of ports in file. At the
//                    end of the list the reference terminal of each group.
//
// t1:g1 t2:g1 ... tn1:g1 t1:g2 t2:g2 .... tnm:gm refg1 refg2 ... refgm

#ifndef NPortT_h
#define NPortT_h 1

class NPortT : public Element
{
	public:

  // The scanner need access to private functions
  friend int NPortT_lex();

  NPortT(const string& iname);

  ~NPortT();

  static const char* getNetlistName()
  {
    return einfo.name;
  }

  // Do some local initialization
  virtual void init() throw(string&);

  // Get a vector with the indexes of the local reference nodes.
  virtual void getLocalRefIdx(UnsignedVector& local_ref_vec,
	TerminalVector& term_list);

  // This element adds equations to the MNAM (in the time domain)
  virtual unsigned getExtraRC(const unsigned& eqn_number,
	const MNAMType& type);
  virtual void getExtraRC(unsigned& first_eqn, unsigned& n_rows) const;

  // Fill MNAM
  virtual void fillMNAM(class FreqMNAM *);
  virtual void fillMNAM(class TimeMNAM *);

	private:

  // Add a new port given the group number
  void addPort(const unsigned& group_number);

  // Set size of t.f. polynomials
  void setNumSize(unsigned tr_idx, unsigned size);
  void setDenSize(unsigned tr_idx, unsigned size);
  // Set polynomials coefficients
  double& setNumCoeff(unsigned tr_idx, unsigned idx);
  double& setDenCoeff(unsigned tr_idx, unsigned idx);

  // Number of ports counted from file
  unsigned ports;
  // Current group number
  unsigned current_group;
  // Vector containing the number of terminals in each group
  UnsignedVector group_vec;

  // Variables used to fill the TimeMNAM
  unsigned ncurr_eqns, nvolt_eqns, my_start_row, extra_rcs;
  DenseIntVector neqns_row, neqns_col;

  // Normalizing frequency
  double fmax;
  // Number of transfer functions
  unsigned ntransfer;
  // vector arrays to hold coefficients
  DenseDoubleVector *numvec, *denvec;
  // Frequency domain workspace matrix
  ComplexDenseMatrix ws_matrix;

  // Element information
  static ItemInfo einfo;

  // Number of parameters of this element
  static const unsigned n_par;

  // Parameter variables
  string filename;
  double mf;

  // Parameter information
  static ParmInfo pinfo[];
};

#endif

