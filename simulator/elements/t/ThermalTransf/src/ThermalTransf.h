// ThermalTransf is a one-port thermal element based on ThermalMMIC
// implementing time transformation.
//
//	          +-----+
//   tref Tamb o--+     +--o 1 pin, tout
//	          +-----+
//
// 1-port description of heatsink mounted MMIC die.
// Default parameters for a L = 400 um^3 GaAs die, with single, central
// surface heating element 0.1L on each side.
//

#ifndef ThermalTransf_h
#define ThermalTransf_h 1

class ThermalTransf : public Element
{
	public:

  ThermalTransf(const string& iname);

  ~ThermalTransf();

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

  // deltaT in time domain
  double getDeltaT(const double& x, const double& ctime, const double& cdt);
  // deltaTheta in transform-domain time (tau)
  int getDeltaTheta(const double& x, const double& tau);

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

  double storedDeltaT, storedDeltaT_new;
  double last_tau, prev_tau, lastTStime, tfrac;
  double Pav, Pold, last_x, bm1;
  int n_tts, n_steps, im1;

  DenseDoubleVector dtheta, storedDeltaT_base, knorm, time;

};

#endif

