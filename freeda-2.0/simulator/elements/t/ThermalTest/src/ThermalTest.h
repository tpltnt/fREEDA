// ThermalTest is a one-port thermal element based on ThermalMMIC
//
//	          +-----+
//   tref Tamb o--+     +--o 1 pin, tout
//	          +-----+
//
// 1-port description of heatsink mounted MMIC die.
// Default parameters for a L = 400 um^3 GaAs die, with single, central
// surface heating element 0.1L on each side.
//

#ifndef ThermalTest_h
#define ThermalTest_h 1

class ThermalTest : public Element
{
	public:

  ThermalTest(const string& iname);

  ~ThermalTest();

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
  double getDeltaT(const double& x);

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
  int Ntimesteps, Nfingers;
  double dt;
  double Tambient;
  bool time_d;
  double L, W, D, xL, xR, yU, yD, Ks, rho, C, b;

  // Internal global variables
  double *Rth, *Pstore;
  double storedDeltaT, storedDeltaT_base, storedDeltaT_old;
  double deltaT_old, deltaT_old1, deltaT_old1_base;
  double currtime, lastTStime, time_old, Pav, Pold;
  int n_ts, n_steps;
};

#endif

