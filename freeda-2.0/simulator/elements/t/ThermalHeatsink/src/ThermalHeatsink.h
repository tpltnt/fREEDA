// This may look like C code, but it is really -*- C++ -*-
//
//
// ThermalHeatsink is an n-port thermal element
//
//	        +-----+
//        0  o--+     |
//        1  o--+     |
//        2  o--+     |
//        3  o--+     |
//        .  o--+     |
//        .  o--+     |
//        .  o--+     |
//           o--+     |
//           o--+     |
//           o--+     |
//        n  o--+     |
//              |     |
// n+1 (ref) o--+     |
//	        +-----+
//
//  N-port description of grid array substrate with NxN, surface
//  heating elements.  Heat loss is entirely by radiation and
//  convection, with no heatsink mounting. Default parameters are for
//  FR-4.  Only a single averaged surface heating element is
//  described, to model the heatspreading effects of heavy
//  metallisation between power dissipating elements.


#ifndef ThermalHeatsink_h
#define ThermalHeatsink_h 1

class ThermalHeatsink : public Element
{
public:

  ThermalHeatsink(const string& iname);

  ~ThermalHeatsink();

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
  double LaplaceDomain(int i,int j, double s);
  double TimeIndependent(int i);
  double_complex LaplaceDomain(int i, double_complex s);
  void TimeDomain(double t, double *Rthtemp);
  double CalculateWeight(int v);
  double fact(int);
  void nrerror(const char*);
  void free_matrix(double** m, int nrl, int nrh, int ncl, int nch);
  double** newmatrix(int nrl, int nrh, int ncl, int nch);
  double *newvector(int, int);
  void free_vector(double*, int, int);

  // Element information
  static ItemInfo einfo;
  // Number of parameters of this element
  static const unsigned n_par;
  // Parameter information
  static ParmInfo pinfo[];
  // Parameter variables
  int Ntimesteps, Ndevices, Narray;
  double dt, Tambient;
  bool time_d, read_input ;
  double L, W, D, xi, eta, epsilon, Ks, rho, C, b, xl, xr, yu, yd;

  // Internal global variables
  int Nsquared;
  double *xL, *xR, *yU, *yD;
  static const double sigma;
  double XL, XR, YU, YD;
  double **Rth, **Pstore, *Rthtemp, *Pstoretotal, *storedDeltaT;

  DenseIntVector my_row;
};

#endif


