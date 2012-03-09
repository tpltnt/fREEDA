// Ymatrix: linear multiport element given by its nodal y parameter matrix.
//          The number of groups and terminals in each group is variable.
//
// Author: original transim version by Mark Summers
//         C++ version and improvements Carlos E. Christoffersen
//
//
//              +-----------+
//         o----|           |
// group1  o----|           |----o  group3
//         o----|           |----o
//         o----|  Ymatrix  |----o
//         |    |           |    |
//         ^    |           |    ^ local reference
//              |           |
// group2  o----|           |----o  ...
//         o----|           |----o  group(m)
//         |    +-----------+    |
//         ^                     ^
//
// Terminal order is: first non-reference terminals of each group. At the
//                    end of the list the reference terminal of each group.
//
// t1:g1 t2:g1 ... tn1:g1 t1:g2 t2:g2 .... tnm:gm refg1 refg2 ... refgm
//
//

#ifndef Ymatrix_h
#define Ymatrix_h 1

class Ymatrix : public Element
{
	public:

  // The scanner need access to private functions
  friend int Ymatrix_lex();

  Ymatrix(const string& iname);

  ~Ymatrix();

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

  // Number of terminals counted from file
  unsigned terms;

  // Current group number
  unsigned current_group;

  // Vector containing the number of terminals in each group
  UnsignedVector group_vec;

  // Current capacity of the internal frequency and matrix vectors
  unsigned vec_capacity;

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
  string pfilename;

  // Parameter information
  static ParmInfo pinfo[];
};

#endif

