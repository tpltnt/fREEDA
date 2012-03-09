// ThermalNportNL is an n+m-port thermal element with non linear boundary conditions
//
//	        +-----+
//        0  o--+     +--o 0
//        1  o--+     +--o 1
//        2  o--+     +--o 2
//        3  o--+     |    .
//        .  o--+     |    .
//        .  o--+     |    .
//        .  o--+     |    .
//           o--+     |    .
//           o--+     |    .
//           o--+     +--o m
//        n  o--+     |
//              |     |
// n+1 (ref) o--+     |
//	        +-----+
//
//  n+m-port description of grid array substrate with NxN surface
//  heating elements and discretised top and bottom surfaces to treat
//  flux non linearity.  Heat loss is entirely by radiation and
//  convection, with no heatsink mounting. Default parameters are for
//  FR-4.  Only a single averaged surface heating element is
//  described, to model the heatspreading effects of heavy
//  metallisation between power dissipating elements.

#ifndef ThermalNportNL_h
#define ThermalNportNL_h 1

class ThermalNportNL : public Element
{
	public:

  ThermalNportNL(const string& iname);

  ~ThermalNportNL();

  static const char* getNetlistName()
  {
    return einfo.name;
  }

  // Do some local initialization
  virtual void init() throw(string&);

  // fill MNAM
  virtual unsigned getExtraRC(const unsigned& eqn_number,
	const MNAMType& type);
  void getExtraRC(unsigned& first_eqn, unsigned& n_rows) const;
  virtual void fillMNAM(FreqMNAM* mnam);

  // State variable transient analysis
  virtual void svTran(TimeDomainSV* tdsv);
  virtual void deriv_svTran(TimeDomainSV* tdsv);

	private:

  // Internal functions
  double LaplaceDomain00(int i, int j, double s);
  double LaplaceDomain0D(int i, int j, double s);
  double LaplaceDomainD0(int i, int j, double s);
  double LaplaceDomainDD(int i, int j, double s);
  double TimeIndependent00(int i, int j);
  double TimeIndependent0D(int i, int j);
  double TimeIndependentD0(int i, int j);
  double TimeIndependentDD(int i, int j);
  double_complex LaplaceDomain00(int i, int j, double_complex s);
  double_complex LaplaceDomain0D(int i, int j, double_complex s);
  double_complex LaplaceDomainD0(int i, int j, double_complex s);
  double_complex LaplaceDomainDD(int i, int j, double_complex s);
  void TimeDomain00(double t, double **Rth00temp);
  void TimeDomain0D(double t, double **Rth0Dtemp);
  void TimeDomainD0(double t, double **RthD0temp);
  void TimeDomainDD(double t, double **RthDDtemp);
  double CalculateWeight(int v);
  double fact(int);
  void nrerror(const char*);
  void free_matrix(double** m, int nrl, int nrh, int ncl, int nch);
  double** newmatrix(int nrl, int nrh, int ncl, int nch);
  double *newvector(int, int);
  void free_vector(double*, int, int);
  double*** newtensor(int nrl, int nrh, int ncl, int nch, int ndl, int ndh);
  void free_tensor(double***, int, int, int, int, int, int);

  // Element information
  static ItemInfo einfo;
  // Number of parameters of this element
  static const unsigned n_par;
  // Parameter information
  static ParmInfo pinfo[];
  // Parameter variables
  int Ntimesteps, Ndevices, Narray, Msubstrate;
  double dt, Tambient;
  bool time_d, read_input ;
  double L, W, D, xi, eta, epsilon, Ks, rho, C, b;

  // Internal global variables
  int Nsquared, Msquared;
  double *xL0, *xR0, *yU0, *yD0, *xLD, *xRD, *yUD, *yDD;
  static const double sigma;
  double ***Rth00, ***Rth0D, ***RthD0, ***RthDD, **Rth00temp, **Rth0Dtemp, **RthD0temp, **RthDDtemp, *Pstoretotal, *storedDeltaT;

  DenseIntVector my_row;
};

#endif

