#include "../../../../analysis/FreqMNAM.h"
#include "../../../../analysis/TimeMNAM.h"
#include "AbmButterbpf10.h"
#include <cstdio>

// This is for a 4-terminal device defined by y = H(s) = I(s) / V(s)
//
//           +------+
//  0 o------+      +------o 2
//     input +      +output
//  1 o------+      +------o 3
//	     +------+
//

// Static members
const unsigned AbmButterbpf10::n_par = 2;

// Element information
ItemInfo AbmButterbpf10::einfo =
{
  "abmbutterbpf10",
  "Two port 10th order Butterworth bandpass filter",
  "Aaron Walker",
  DEFAULT_ADDRESS"category:behavioral",
  "2003_05_15"
};

// Parameter information
ParmInfo AbmButterbpf10::pinfo[] =
{
  {"bw", "Filter 3dB bandwidth", TR_DOUBLE, true},
  {"f0", "Filter center frequency", TR_DOUBLE, true}
};


AbmButterbpf10::AbmButterbpf10(const string& iname) : Element(&einfo, pinfo, n_par, iname)
{
  // Set default parameter values
  paramvalue[0] = &(bw = 1);
  paramvalue[1] = &(f0 = 1);

  // Set the number of terminals
  setNumTerms(4);
  // Set flags
  setFlags(LINEAR | ONE_REF | TR_FREQ_DOMAIN);
}

// init() function
void AbmButterbpf10::init() throw(string&)
{
  unsigned int i;
  // Initialize coefficient vectors
  k0.resize(5);
  k1.resize(5);
  bt.resize(5);
  a.resize(5);
  b.resize(5);
  c.resize(5);
  d.resize(5);
  e.resize(5);
  f.resize(5);
  h.resize(5);

  // Rational fraction residues for lowpass prototype
  k0[0] = 0.;
  k0[1] = 19.5134116599;
  k0[2] = -0.37174803446;
  k0[3] = 3.7977275657;
  k0[4] = -22.93939119;

  k1[0] = -8.690133435;
  k1[1] = 23.474671498;
  k1[2] = 0.4472135955;
  k1[3] = 2.823595516;
  k1[4] = -17.05534717;

  // Rational fraction first order coefficients from lowpass
  // prototype.  Second and zero order coefficients are 1.
  bt[0] = 1.41421356;
  bt[1] = 1.97537668;
  bt[2] = 0.31286893;
  bt[3] = 0.907981;
  bt[4] = 1.78201305;

  // Coefficients for numerator and denominator polynomials
  // after bandpass transform.  Each rational fraction has
  // the following form.
  //
  //    H_bp(s) =        a*s^2 + b*s + c
  //              -----------------------------
  //              s^4 + d*s^3 + e*s^2 + f*s + h
  //

  // Some of the coefficients are constant for all of the
  // fractions while others depend on bt, k0, and k1.

  bw_rad  = bw * 2. * pi;
  w0  = f0 * 2. * pi;
  w02 = w0 * w0;
  w04 = w02 * w02;

  cout << endl;
  for(i=0;i < bt.length(); ++i)
	{
    a[i] = k0[i]*bw_rad;
    b[i] = k1[i]*bw_rad*bw_rad;
    c[i] = k0[i]*bw_rad*w02;
    d[i] = bt[i]*bw_rad;
    f[i] = bt[i]*bw_rad*w02;
    printf("a%d= %f, b%d= %f, c%d= %f, d%d= %f, f%d= %f\n",i,a[i],i,b[i],i,c[i],i,d[i],i,f[i]);
  }
  e[1] = 2*w02+bw_rad*bw_rad;
  h[1] = w04;

  test.resize(5);
  test[0] = 5;
  test[1] = 3;
}

unsigned AbmButterbpf10::getExtraRC(const unsigned& eqn_number, const MNAMType& type)
{
  unsigned extra_rcs;
  if (type == TIME_DOMAIN)
	{
    // Add extra rows and columns, the extra rows and columns necessary
    // correspond to the derivatives of the time domain current and all
    // derivatives above order one in the time domain voltage.  These
    // terms are necessary from the transform of the general transfer
    // function in the s domain.
    //
    // Terms needed in u vector of nodal voltages:
    //
    //
    //   u(t) = |v0|
    //          |v1|
    //          |...|
    //		|i|
    //		|di/dt|
    //		|d^2i/dt^2|
    //		|d^3i/dt^3|
    // 		|dv0/dt|
    //          |dv1/dt|
    //          |d^2v0/dt^2|
    // 		|d^2v1/dt^2|

    my_start_row = eqn_number;
    // extra  = # of current terms + 2*(# of voltage terms added)
    // there are four current terms added for each rational fraction
    // zero-third derivatives and four voltage terms for each rational
    // fraction, first and second derivatives for two terminals.
    extra_rcs = 4*bt.length() + 2*(b.length()+ c.length());
		//    cout << endl << "extra_rcs " << extra_rcs << endl;
		//    cout << endl << "my_start_row " << my_start_row << endl;
  }
  else
	{
    // No need for extra rows and columns in frequency domain
    my_start_row = 0;
    extra_rcs = 0;
  }
  return extra_rcs;
}

void AbmButterbpf10::getExtraRC(unsigned& first_eqn, unsigned& n_rows) const
{
  first_eqn = my_start_row;
  n_rows = extra_rcs;
}

void AbmButterbpf10::fillMNAM(FreqMNAM* mnam)
{
  double_complex jw = mnam->getFreq() * double_complex(0., 2.) * pi;
  double_complex jw2, jw3, jw4;
  unsigned int i;

  // Evaluate H(jw)
  double_complex g = 0;
  cout << endl << "Freq MNAM"<< endl;
  jw2 = jw*jw;
  jw3 = jw*jw2;
  jw4 = jw*jw3;
  for(i = 0; i < a.length(); ++i)
  {
    g += (jw3*a[i] + jw2*b[i]+ jw*c[i]) /
		(jw4 + jw3*d[i] + jw2*e[1] + jw*f[i] + h[1]);
  }

  // Set the frequency domain stamp matrix.  Rows are currents here,
  // currents at terminals 2 and 3.  Columns are voltages, here only
  // considering forward transfer functions so columns are 0 and 1.
  if(g != double_complex(0))
	{
    mnam->setQuad(getTerminal(2)->getRC(), getTerminal(3)->getRC(),
		getTerminal(0)->getRC(), getTerminal(1)->getRC(),g);
  }
}


void AbmButterbpf10::fillMNAM(TimeMNAM* mnam)
{
  unsigned int i;
	int j;
  int extra_rows_each = 8;
  // Additional rows for M (G) matrix
  cout << endl << "Time MNAM"<< endl;

  for (i = 0;i < b.length(); i++)
	{
    // Add column terms to relate added current term to i2 and i3
    // iadd = -i2 = i3
    mnam->setMElement(getTerminal(2)->getRC(), my_start_row + i*extra_rows_each, -1.);
    mnam->setMElement(getTerminal(3)->getRC(), my_start_row + i*extra_rows_each, 1.);

    // Add term for coefficient h for each rational fraction in the transfer
    // function.   Gu = sv_part_1
    mnam->setMElement(my_start_row + i*extra_rows_each,
		my_start_row + i*extra_rows_each, h[1]);

    // Add terms to include extra elements of the nodal voltage vector u
    for(j=1;j < extra_rows_each; j++){
      mnam->setMElement(my_start_row + i*extra_rows_each + j,
			my_start_row + i*extra_rows_each + j, 1.);
    }
  }

  // Additional rows for Mp (C) matrix

  // Add row for constitutive relation for transfer function involving
  // du/dt terms. Cdu/dt = sv_part_2
  for (i = 0; i < b.length(); i++)
	{
    // Add term for di/dt term in du(t)/dt
    mnam->setMpElement(my_start_row + i*extra_rows_each,
		my_start_row + i*extra_rows_each, f[i]);
    // Add term for di^2/dt^2 term
    mnam->setMpElement(my_start_row + i*extra_rows_each,
		my_start_row + i*extra_rows_each + 1, e[1]);
    // Add term for di^3/dt^3 term
    mnam->setMpElement(my_start_row + i*extra_rows_each,
		my_start_row + i*extra_rows_each + 2, d[i]);
    // Add term for di^4/dt^4 term
    mnam->setMpElement(my_start_row + i*extra_rows_each,
		my_start_row + i*extra_rows_each + 3, 1.);
    // Add term for dv0/dt term
    mnam->setMpElement(my_start_row + i*extra_rows_each,
		getTerminal(0)->getRC(), -c[i]);
    // Add term for dv1/dt term
    mnam->setMpElement(my_start_row + i*extra_rows_each,
		getTerminal(1)->getRC(), c[i]);
    // Add term for dv0^2/dt^2 term
    mnam->setMpElement(my_start_row + i*extra_rows_each,
		my_start_row + i*extra_rows_each + 4, -b[i]);
    // Add term for dv1^2/dt^2 term
    mnam->setMpElement(my_start_row + i*extra_rows_each,
		my_start_row + i*extra_rows_each + 5, b[i]);
    // Add term for dv0^3/dt^3 term
    mnam->setMpElement(my_start_row + i*extra_rows_each,
		my_start_row + i*extra_rows_each + 6, -a[i]);
    // Add term for dv1^3/dt^3 term
    mnam->setMpElement(my_start_row + i*extra_rows_each,
		my_start_row + i*extra_rows_each + 7, a[i]);

    // Equate derivative terms in du/dt to those in u,
    // i.e. di/dt(u) - di/dt(du/dt) = 0

    // Equate di/dt
    mnam->setMpElement(my_start_row + i*extra_rows_each + 1,
		my_start_row + i*extra_rows_each, -1.);
    // Equate di^2/dt^2
    mnam->setMpElement(my_start_row + i*extra_rows_each + 2,
		my_start_row + i*extra_rows_each + 1, -1.);
    // Equate di^3/dt^3
    mnam->setMpElement(my_start_row + i*extra_rows_each + 3,
		my_start_row + i*extra_rows_each + 2, -1.);
    // Equate dv0/dt
    mnam->setMpElement(my_start_row + i*extra_rows_each + 4,
		getTerminal(0)->getRC(), -1.);
    // Equate dv1/dt
    mnam->setMpElement(my_start_row + i*extra_rows_each + 5,
		getTerminal(1)->getRC(), -1.);
    // Equate dv0^2/dt^2
    mnam->setMpElement(my_start_row + i*extra_rows_each + 6,
		my_start_row + i*extra_rows_each + 4, -1.);
    // Equate dv1^2/dt^2
    mnam->setMpElement(my_start_row + i*extra_rows_each + 7,
		my_start_row + i*extra_rows_each + 5, -1.);
  }
}

