// ThermalHeatsink1 is a one-port thermal element
//
//	          +-----+
//   tref Tamb o--+     +--o 1 pin, tout
//	          +-----+
//
//  1-port description of grid array substrate with a single,
//  averaged, surface heating element.  Heat loss is entirely by
//  radiation and convection, with no heatsink mounting.  Default
//  parameters are for FR-4.  Only a single averaged surface heating
//  element is described, to model the heatspreading effects of heavy
//  metallisation between power dissipating elements.
//

#ifndef ThermalHeatsink1_h
#define ThermalHeatsink1_h 1

class ThermalHeatsink1 : public Element
{
	public:

  ThermalHeatsink1(const string& iname);

  ~ThermalHeatsink1();

  static const char* getNetlistName()
  {
    return einfo.name;
  }

  // Do some local initialization
  virtual void init() throw(string&);

  // fill MNAM
  virtual void fillMNAM(FreqMNAM* mnam);

  // State variable transient analysis
  virtual void svTran(TimeDomainSV* tdsv);
  virtual void deriv_svTran(TimeDomainSV* tdsv);

	private:

  // Internal functions
  double LaplaceDomain(double Tambient, double s);
  double_complex CLaplaceDomain(double Tambient, double_complex s);
  double TimeDomain(double Tambient, double t);
  double TimeIndependent(double Tambient);
  double CalculateWeight(int v);
  double fact(int);
  void nrerror(const char*);
  double *newvector(int, int);
  void free_vector(double*, int, int);

  // Element information
  static ItemInfo einfo;
  // Number of parameters of this element
  static const unsigned n_par;
  // Parameter information
  static ParmInfo pinfo[];
  // Parameter variables
  int Ntimesteps, Ndevices;
  double dt;
  double Tambient;
  bool time_d;
  double L, W, D, xL, xR, yU, yD, Ks, rho, C, xi, eta, epsilon, b;

  // Internal global variables
  double *Rth, *Pstore;
  double storedDeltaT;

};

#endif

