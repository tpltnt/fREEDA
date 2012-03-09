// ThermalHeatsinkMMIC1 is a one-port thermal element
//
//	          +-----+
//   tref Tamb o--+     +--o 1 pin, tout
//	          +-----+
//
// 1-port description of heatsink mounted MMIC die.
// Default parameters for a L = 400 um^3 GaAs die, with single, central
// surface heating element 0.1L on each side.
//

#ifndef ThermalHeatsinkMMIC1_h
#define ThermalHeatsinkMMIC1_h 1

class ThermalHeatsinkMMIC1 : public Element
{
	public:

  ThermalHeatsinkMMIC1(const string& iname);

  ~ThermalHeatsinkMMIC1();

  static const char* getNetlistName()
  {
    return einfo.name;
  }

  // Do some local initialization
  virtual void init() throw(string&);

  // Functions used for zero SV approach
  unsigned getExtraRC(const unsigned& eqn_number, const MNAMType& type);
  void getExtraRC(unsigned& first_eqn, unsigned& n_rows) const;
  void fillMNAM(TimeMNAM* mnam);
  void setLastResult(DenseDoubleVector& res, const double& time);
  void fillSourceV(TimeMNAM* mnam);

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
  int Ntimesteps, Nfingers;
  double dt;
  double Tambient;
  bool time_d, zero_sv, kt;
  double L, W, D, xL, xR, yU, yD, Ks, rho, C, b;

  // Internal global variables
  double *Rth, *Pstore;
  double storedDeltaT, storedTime;
  int nstep;
  unsigned my_row;
};

#endif

