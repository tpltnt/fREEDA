// NPort: linear multiport element given by its y parameter matrix.
//        The number of groups and terminals in each group is variable.
//
// Author: Carlos E. Christoffersen
//         Based on original code by Mark Summers.
//
//
//              +-----------+
//         o----|           |
// group1  o----|           |----o  group3
//         o----|           |----o
//         o----|  NPort    |----o
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
//
//

#ifndef NPort_h
#define NPort_h 1

class NPort : public Element
{
	public:

  // The scanner need access to private functions
  friend int NPort_lex();

  NPort(const string& iname);

  ~NPort();

  static const char* getNetlistName()
  {
    return einfo.name;
  }

  // Check is element has been built correctly
  virtual void check() throw(string&);

  // Do some local initialization
  virtual void init() throw(string&);


  // Get a vector with the indexes of the local reference nodes.
  virtual void getLocalRefIdx(UnsignedVector& local_ref_vec,
	TerminalVector& term_list);

  // Fill MNAM
  virtual void fillMNAM(class FreqMNAM *);

	private:

  // Add a new port given the group number
  void addPort(const unsigned& group_number);

  // Add one frequency and matrix
  ComplexDenseMatrix* newMatrix(const double& frequency);

  // Workspace matrix
  ComplexDenseMatrix ws_matrix;

  // Number of ports counted from file
  unsigned ports;

  // Current group number
  unsigned current_group;

  // Vector containing the number of terminals in each group
  UnsignedVector group_vec;

  // Current size of the internal frequency and matrix vectors
  unsigned vec_size;

  // Internal frequency vector
  double* freq_vec;

  // Internal y paramenter matrix vector
  ComplexDenseMatrix** ymat_vec;

  // Last frequency index.
  int last_index;

  // Element information
  static ItemInfo einfo;

  // Number of parameters of this element
  static const unsigned n_par;

  // Parameter variables
  string filename;
  int max_freq;

  // Parameter information
  static ParmInfo pinfo[];

};

#endif

