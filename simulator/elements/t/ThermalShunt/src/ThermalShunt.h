// This may look like C code, but it is really -*- C++ -*-
//
//
// ThermalShunt is a two-port thermal element
//
//	           +-----+     
//     t0     0 o--+     |
//                 |     |
//     tD     1 o--+     |
//                 |     |
//     Tref   2 o--+     |
//	           +-----+     
//
// 2-port description of MMIC die with variable base temperature,
// vertical matching at interfaces.  
// parameters for a L = 400 um^3 GaAs die, with single, central
// surface heating element 0.1L on each side. 
//

#ifndef ThermalShunt_h
#define ThermalShunt_h 1

class ThermalShunt : public Element
{
public:
  
  ThermalShunt(const string& iname);
  
  ~ThermalShunt();

  static const char* getNetlistName()
  {
    return einfo.name;
  }

  // Do some local initialization
  virtual void init() throw(string&);

  // This element adds equations to the MNAM
  virtual unsigned getExtraRC(const unsigned& eqn_number, 
			      const MNAMType& type);
  virtual void getExtraRC(unsigned& first_eqn, unsigned& n_rows) const;
  // fill MNAM
  virtual void fillMNAM(FreqMNAM* mnam);

  // State variable transient analysis
  virtual void svTran(TimeDomainSV* tdsv);
  virtual void deriv_svTran(TimeDomainSV* tdsv);

private:

  // Private function list
  double LaplaceDomain00(double s);
  double LaplaceDomain0D(double s);
  double LaplaceDomainD0(double s);
  double LaplaceDomainDD(double s);
  double_complex CLaplaceDomain00(double_complex s);
  double_complex CLaplaceDomain0D(double_complex s);
  double_complex CLaplaceDomainD0(double_complex s);
  double_complex CLaplaceDomainDD(double_complex s);
  double TimeIndependent00();
  double TimeDomain00(double t);
  double TimeDomain0D(double t);
  double TimeDomainD0(double t);
  double TimeDomainDD(double t);
  // I got all there unchanged from Thermal1.cc
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
  // Parameter and internal variables
  int Ntimesteps, Nfingers;
  double dt, L, W, D, xL, xR, yU, yD, Pnorm0, PnormD;
  double xi, eta, sigma, epsilon, T0, Tambient, H, Ks, rho, C, k, b;
  bool time_d, read_input;

  double *Rth00, *Rth0D, *RthD0, *RthDD, *P0store, *PDstore;
  double storedDeltaT0, storedDeltaTD;

  int my_row1, my_row2;
};

#endif


