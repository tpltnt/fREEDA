#include "../../../../network/ElementManager.h"
#include "../../../../analysis/TimeDomainSV.h"
#include <cmath>
#include "Zchebyi.h"
extern "C"
{
#include "../../../../inout/ftvec.h"
#include "../../../../inout/report.h"
}

// Static members
const unsigned Zchebyi::n_par = 11 ;

// Element information
ItemInfo Zchebyi::einfo =
{
  "zchebyi",
  "Z-Domain Chebyshev Type I Lowpass or Bandpass Filter",
  "Frank P. Hart",
  ""
};

// Parameter information
ParmInfo Zchebyi::pinfo[] =
{
  // Netlist must specify flp&fhp OR fo&pbw, but not both
  {"flp", "Lower Pass Frequency (Hz)",TR_DOUBLE, false},
  {"fhp", "Upper Pass Frequency (Hz)",TR_DOUBLE, false},
  {"fo", "Center Frequency (Hz)",TR_DOUBLE, false},
  {"pbw", "Passband bandwidth (Hz or \% if less than 100)",TR_DOUBLE, false},
  // pbfdb is optional, -3 assumed if not set
  {"pbfdb", "Passband Flatness (negative) (dB)",TR_DOUBLE, false},
  // Netlist must specify fls&fhs&sbadb OR pord, but not both
  {"fls", "Lower Stop Frequency (Hz)",TR_DOUBLE, false},
  {"fhs", "Upper Stop Frequency (Hz)",TR_DOUBLE, false},
  {"sbadb", "Stopband Attenuation (negative) (dB)",TR_DOUBLE, false},
  {"pord", "Lowpass Prototype Filter Order",TR_INT, false},
  {"ildb", " Insertion Loss (negative) (dB)",TR_DOUBLE, false},
  {"rep", "Record filter design info in report()",TR_BOOLEAN, false},
};


Zchebyi::Zchebyi(const string& iname) : Element(&einfo, pinfo, n_par, iname)
{
  // Value required for fhp and fhs in netlist; defaults for others
  paramvalue[0] = &(flp=0.0) ;
  paramvalue[1] = &(fhp=0.0) ;
  paramvalue[2] = &(fo=0.0) ;
  paramvalue[3] = &(pbw=0.0) ;
  paramvalue[4] = &(pbfdb=-3.0) ;
  paramvalue[5] = &(fls=0.0) ;
  paramvalue[6] = &(fhs=0.0) ;
  paramvalue[7] = &(sbadb=-10.0) ;
  paramvalue[8] = &(pord=0);
  paramvalue[9] = &(ildb=0.0);
  paramvalue[10] = &(rep=false) ;

  setFlags(NONLINEAR | MULTI_REF | TR_TIME_DOMAIN); // Set flags
  setNumTerms(4);	// Set the number of terminals
  setNumberOfStates(2);	// Set number of states

  name_len = iname.length();
  name_ch = new char[name_len + 1];
  iname.copy(name_ch, name_len, 0);
  name_ch[name_len] = 0;
}

Zchebyi::~Zchebyi() // 'delete' in reverse order of 'new'
{
  if (lowpass)
  {
    delete [] y_int;
    delete [] x_int;
    delete [] w;
    delete [] w_delp;
    delete [] w_del2;
    delete [] w_del1;
    delete [] z_den3;
    delete [] z_den2;
    delete [] z_den1;
    delete [] z_num3;
    delete [] z_num2;
    delete [] z_num1;
  }
  else
  { // bandpass...
    delete [] y_int;
    delete [] x_int;
    delete [] w;
    delete [] w_delp;
    delete [] w_del4;
    delete [] w_del3;
    delete [] w_del2;
    delete [] w_del1;
    delete [] z_den5;
    delete [] z_den4;
    delete [] z_den3;
    delete [] z_den2;
    delete [] z_den1;
    delete [] z_num3;
    delete [] z_num2;
    delete [] z_num1;
  }
  delete [] plp_den3;
  delete [] plp_den2;
  delete [] plp_den1;
  delete [] name_ch;
}

void Zchebyi::getLocalRefIdx(UnsignedVector& local_ref_vec,
    TerminalVector& term_list)
{
  // Make sure the vectors are empty
  term_list.erase(term_list.begin(), term_list.end());
  local_ref_vec.erase(local_ref_vec.begin(), local_ref_vec.end());

  // Insert vector elements
  term_list.push_back(getTerminal(0)); // Input port terminal
  term_list.push_back(getTerminal(1)); // Input port local reference terminal
  local_ref_vec.push_back(1); // Local reference index

  term_list.push_back(getTerminal(2)); // Output port terminal
  term_list.push_back(getTerminal(3)); // Output port local reference terminal
  local_ref_vec.push_back(3); // Local reference index
}

void Zchebyi::init() throw(string&)
{
  if(isSet(&pbfdb) && (pbfdb>0.0))
  {
    throw(" Passband Flatness must be negative! ") ;
  }
  else if(isSet(&sbadb) && (sbadb>0.0))
  {
    throw(" Stopband Attenuation must be negative! ") ;
  }
  else if(isSet(&ildb) && (ildb>0.0))
  {
    throw(" Insertion Loss must be negative! ") ;
  }
  //else if( (isSet(&flp)||isSet(&fhp)) && (isSet(&fo)||isSet(&pbw)) ){
  //  throw(" isSet detected flp,fhp // fo,pbw conflict! ") ;
  //}
  //else if( (isSet(&fls)||isSet(&fhs)) && isSet(&pord) ) {
  //  throw(" isSet detected fls,fhs // pord conflict! ") ;
  //}
  //else if( ((flp!=0.0)||(fhp!=0.0)) && ((fo!=0.0)||(pbw!=0.0)) ){
  //  throw(" Specify flp & fhp OR fo & pbw, but not both! ") ;
  //}
  //else if( ((fls!=0.0)||(fhs!=0.0)) && (pord!=0)) {
  //  throw(" Specify fls & fhs OR pord, but not both! ") ;
  //}
  else
  {
    Ngiven = false;
    if (pord!=0) Ngiven = true;
    do_filter_init = true; // init filter coeffs/delays in first eval
    x_last = zero;	// change to 0?
    y_last = zero;	// change to 0?
    x[0] = zero; x[1] = zero;
    vp[0] = zero; vp[1] = zero;
    ip[0] = zero; ip[1] = zero;
    //name_str = &iname;
    //name_len = iname.length();
  }
}

// This routine is called at each time step
// (Actually, several times at each time step because the Jacobian
// is evaluated numerically).
void Zchebyi::svTran(TimeDomainSV* tdsv)
{
  if(tdsv->DC())
  {
    tdsv->u(0) = tdsv->getX(0);
    tdsv->i(0) = zero;
    tdsv->i(1) = tdsv->getX(1);
    tdsv->u(1) = zero;
  }
  else
  {
    first_iteration = tdsv->firstTime();
    if (first_iteration) { // new first iter, so *can't* set x_last to
      x_last = zero;      // x[0] because the analysis routine uses
      y_last = zero;      // same value as new time step 1st iterate!
    }
    else
    {
      x_last = x[0]; // save previous iterate values before
      y_last = vp[1];
    }
    // The state vars are (0) the in. volt.; (1) out. (free) current
    // Note: If can eliminate the free state var, so much the better..
    tdsv->u(0) = tdsv->getX(0);
    tdsv->i(1) = tdsv->getX(1);
    x[0] = tdsv->u(0);
    x[1] = tdsv->i(1);

    eval(x, vp, ip, tdsv);

    // Set output variables...
    tdsv->u(1) = vp[1];
    tdsv->i(0) = ip[0];
  }
}

void Zchebyi::deriv_svTran(TimeDomainSV* tdsv)
{
  tdsv->getJu()(0,0) = one;
  tdsv->getJu()(0,1) = zero;
  //tdsv->getJu()(1,0) = zero;
  tdsv->getJu()(1,1) = one;

  tdsv->getJi()(0,0) = zero;
  tdsv->getJi()(0,1) = zero;
  tdsv->getJi()(1,0) = zero;
  tdsv->getJi()(1,1) = one;
  if(tdsv->DC())
  {
    if ( ((vp[1]-y_last)<=zero) && ((x[0]-x_last)>=zero) )
      // if ( ((vp[1]-y_last)<=zero) && ((x[0]-x_last)<=zero) )
      tdsv->getJu()(1,0) = zero;
    else
      tdsv->getJu()(1,0) = (vp[1] - y_last)/(x[0] - x_last);
  }
  else
  {
    assert(!tdsv->firstTime());
    if ( ((vp[1]-y_last)<=zero) && ((x[0]-x_last)>=zero) )
      // if ( ((vp[1]-y_last)<=zero) && ((x[0]-x_last)<=zero) )
      tdsv->getJu()(1,0) = zero;
    else
      tdsv->getJu()(1,0) = (vp[1] - y_last)/(x[0] - x_last);
    // jacval = tdsv->getJu()(1,0);
  }
}

void Zchebyi::eval(double* x, double* vp,
    double* ip, const TimeDomainSV* tdsv)
{
  // init. filter on first eval entry (so timestep is valid)
  if (do_filter_init)
  {
    tee = tdsv->getdt();  // set tee to valid number...
    filter_init();        // figure out filter form,coeffs,delays...
    do_filter_init = false; // don't do this a second time
  }

  // Are we at first iteration in time step?
  first_iteration = tdsv->firstTime();

  // State variable assignments to element physical quantities...
  vp[0] = x[0] ;	// x[0] = v_in(t)
  ip[1] = x[1] ;        // x[1] = i_out(t)

  // Scale by constant to make DC response 0, then proceed...
  x_int[kmax-1] = z_const*x[0];

  if (lowpass)
  {
    if (first_iteration)
    {	// just ended newton it, update delays
      for (k=kmax-1;k>=0;k--)
      {
        w_del2[k] = w_del1[k];
        w_del1[k] = w_delp[k];
      }
    }
    for (k=kmax-1;k>=0;k--)
    {
      w[k] = x_int[k] - ( z_den2[k]*w_del1[k] + z_den3[k]*w_del2[k] );
      y_int[k] =  z_num1[k]*w[k] + z_num2[k]*w_del1[k]
        + z_num3[k]*w_del2[k];
      w_delp[k] = w[k];  // save present iterate value
      if (k>0)
      {
        x_int[k-1] = y_int[k];  // pass along to next block
      }
    }
  }
  else { // bandpass
    if (first_iteration)
    {	// just ended newton it, update delays
      for (k=kmax-1;k>=0;k--)
      {
        w_del4[k] = w_del3[k];
        w_del3[k] = w_del2[k];
        w_del2[k] = w_del1[k];
        w_del1[k] = w_delp[k];
      }
    }
    for (k=kmax-1;k>=0;k--)
    {
      w[k] = x_int[k] - ( z_den2[k]*w_del1[k] + z_den3[k]*w_del2[k]
          + z_den4[k]*w_del3[k] + z_den5[k]*w_del4[k] );
      y_int[k] =  z_num1[k]*w[k] + z_num2[k]*w_del2[k]
        + z_num3[k]*w_del4[k];
      w_delp[k] = w[k];   // save present iterate value
      if (k>0)
      {
        x_int[k-1] = y_int[k];  // pass along to next block
      }
    }
  }

  //Assign filter output value to vp[1]...
  ip[0] = zero ;      // Should this be 0?  It's not really 0...
  vp[1] = y_int[0] ;  // Use when y_int is adouble.

}

double Zchebyi::acosh(double z)
{
  return log(z + pow((z*z - 1.0),0.5));
}

double Zchebyi::asinh(double z)
{
  return log(z + pow((z*z + 1.0),0.5));
}

void Zchebyi::filter_init()
{
  // end of fREEDA-typical stuff, rest is Z-domain internals setup

  lowpass = false; 	// assume bandpass by default
  if ((fo==0.0)&&(pbw==0.0)&&(flp==0.0))
  {
    lowpass = true;
  }

  if (~lowpass)
  { // figure out passband info...
    if ((fo!=0.0)&(pbw!=0.0))
    {
      if (pbw<100)
      { // assume pbw is passband percentage of fo
        pbw = (pbw/100)*fo;
      }
      flp = -pbw/2.0 + sqrt((pbw*pbw/4.0) + fo*fo);
      fhp = pbw/2.0 + sqrt((pbw*pbw/4.0) + fo*fo);
    }
  }
  if (Ngiven)
  {
    sbadb=-40;
    N = pord;
  }

  // convert specs in dB to absolutes...
  if (pbfdb==-3.0)
  {
    epsilon = 1.0;
  }
  else
  {
    epsilon = sqrt(pow(10.0,-pbfdb/10.0)-1.0);
  }
  lambda = sqrt(pow(10.0,-sbadb/10.0)-1.0);

  if (Ngiven)
  { // must determine sbw for Cheby filters... 1/24/05
    sbw = pbw*pow(10.0,1/N)*(lambda/epsilon);
    fls = -sbw/2.0 + sqrt((sbw*sbw/4.0) + fo*fo);
    fhs = sbw/2.0 + sqrt((sbw*sbw/4.0) + fo*fo);
  }
  else
  {
    sbw = 0.0;
  }

  // convert to radian frequencies and prewarp...
  wls = tan(pi*fls*tee);
  wlp = tan(pi*flp*tee);
  whs = tan(pi*fhs*tee);
  whp = tan(pi*fhp*tee);

  // determine prototype LP filter's critical freqs...
  if (lowpass)
  {
    wp = whp;
    ws = whs;
    wo2 = 0.0;
    wo = 0.0;
  }
  else
  {
    wo2 = wlp*whp;
    wo = sqrt(wo2);
    wp = (whp*whp - wo2)/whp;
    wsh = (whs*whs - wo2)/whs;
    if (wls == 0.0)
    {
      wsl = wsh;
    }
    else
    {
      wsl = (wo2 - wls*wls)/wls;
    }
    if (wsh > wsl)
    {
      ws = wsh;
    }
    else
    {
      ws = wsl;
    }
  }

  // determine proto LP's order, N, and cutoff freq...
  if (Ngiven==false)
  {
    N = static_cast<int>(ceil(acosh(lambda/epsilon)/acosh(ws/wp)));
  }
  wc = wp;

  // determine whether filter order even or odd
  // Even order -- pole pairs only
  // Odd order -- pole pairs plus single pole at (normalized) -1
  if (N%2 == 0)
  {
    Nodd = false;
  }
  else
  {
    Nodd = true;
  }

  // determine maximum k for filter coeff&delay 2-D arrays
  if (Nodd)
  {
    kmax = (N+1)/2;
  }
  else
  {
    kmax = N/2;
  }

  // 1/24/05
  double vo = asinh(1.0/epsilon)/N;

  // initialize storage for proto LP filter coeffs, etc...
  plp_den1 = new double[kmax];
  plp_den2 = new double[kmax];
  plp_den3 = new double[kmax];
  for (k=0;k<kmax;k++)
  {
    plp_den1[k] = 1.0;
    plp_den2[k] = wc*2.0*sinh(vo)*sin((2.0*k+1)*pi/(2.0*N));
    plp_den3[k] = wc*wc*
      (sinh(vo)*sinh(vo)*
       sin((2.0*k+1)*pi/(2.0*N))*sin((2.0*k+1)*pi/(2.0*N))
       +cosh(vo)*cosh(vo)*
       cos((2.0*k+1)*pi/(2.0*N))*cos((2.0*k+1)*pi/(2.0*N))
      );
  }
  if (Nodd)
  {
    plp_den1[kmax-1] = 0.0;
    plp_den2[kmax-1] = 1.0;
    plp_den3[kmax-1] = wc*sinh(vo);
  }

  // Set proto LP filter for DC gain of 1
  // 1/24/05
  if (Nodd)
    plp_num = -1.0;
  else
    plp_num = 1.0/sqrt(1.0 + epsilon*epsilon);
  for (k=0;k<kmax;k++)
    plp_num = plp_num*plp_den3[k];

  // Include insertion loss effects...
  if (ildb != 0.0)
  {
    iltf = sqrt(pow(10.0,ildb/10.0));
  }
  else
  {
    iltf = 1.0;
  }

  // Initialize Z constant mult to proto LP result...
  z_const = plp_num*iltf;

  // Start the z-filter definition coding...
  if (lowpass)
  {
    // initialize storage for Z numerator coeffs, etc...
    z_num1 = new double[kmax]; // coefficient of z^0
    z_num2 = new double[kmax]; // coefficient of z^1
    z_num3 = new double[kmax]; // coefficient of z^2
    for (k=0;k<kmax;k++)
    {
      z_num1[k] = 1.0;
      z_num2[k] = 2.0;
      z_num3[k] = 1.0;
    }
    if (Nodd)
    {
      z_num1[kmax-1] = 1.0;
      z_num2[kmax-1] = 1.0;
      z_num3[kmax-1] = 0.0;
    }

    // initialize storage for Z denominator coeffs, etc...
    z_den1 = new double[kmax]; // coefficient of z^0
    z_den2 = new double[kmax]; // coefficient of z^1
    z_den3 = new double[kmax]; // coefficient of z^2
    for (k=0;k<kmax;k++)
    {
      a = plp_den1[k]; b = plp_den2[k]; c = plp_den3[k];
      z_den1[k] = a + b + c;
      z_den2[k] = - 2.0*a + 2.0*c;
      z_den3[k] = a - b + c;
    }
    if (Nodd)
    {
      z_den1[kmax-1] = b + c;
      z_den2[kmax-1] = -b + c;
      z_den3[kmax-1] = 0.0;
    }

    // Normalize denominator coeffs to ease filter implementation
    for (k=0;k<kmax;k++)
    {
      if ((z_den1[k] != 0.0) && (z_den1[k] != 1.0))
      {
        z_const = z_const/z_den1[k];
        z_den2[k] = z_den2[k]/z_den1[k];
        z_den3[k] = z_den3[k]/z_den1[k];
        z_den1[k] = 1.0;
      }
    }

    // initialize storage for delay elements, etc...
    w_del1 = new double[kmax]; // first delay element
    w_del2 = new double[kmax]; // second delay element
    w_delp = new double[kmax]; // present element (for iteration)
    for (k=0;k<kmax;k++)
    {
      w_del1[k] = 0.0;
      w_del2[k] = 0.0;
      w_delp[k] = 0.0;
    }

    // initialize intermediate storage for debug tracking...
    w = new double[kmax];
    x_int = new double[kmax];
    y_int = new double[kmax];
    for (k=0;k<kmax;k++)
    {
      w[k] = 0.0;
      x_int[k] = 0.0;
      y_int[k] = 0.0;
    }
  }
  else
  { // bandpass
    // initialize storage for Z numerator coeffs, etc...
    z_num1 = new double[kmax]; // coefficient of z^0
    z_num2 = new double[kmax]; // coefficient of z^2
    z_num3 = new double[kmax]; // coefficient of z^4
    for (k=0;k<kmax;k++)
    {
      z_num1[k] = 1.0;
      z_num2[k] = -2.0;
      z_num3[k] = 1.0;
    }
    if (Nodd)
    {
      z_num1[kmax-1] = 1.0;
      z_num2[kmax-1] = -1.0;
      z_num3[kmax-1] = 0.0;
    }

    // initialize storage for Z denominator coeffs, etc...
    z_den1 = new double[kmax]; // coefficient of z^0
    z_den2 = new double[kmax]; // coefficient of z^1
    z_den3 = new double[kmax]; // coefficient of z^2
    z_den4 = new double[kmax]; // coefficient of z^3
    z_den5 = new double[kmax]; // coefficient of z^4
    for (k=0;k<kmax;k++)
    {
      a = plp_den1[k]; b = plp_den2[k]; c = plp_den3[k];
      z_den1[k] = a + b + 2.0*wo2*a + c + b*wo2 + a*wo2*wo2;
      z_den2[k] = - 4.0*a - 2.0*b + 2.0*b*wo2 + 4.0*a*wo2*wo2;
      z_den3[k] = 6.0*a - 4.0*wo2*a - 2.0*c + 6.0*a*wo2*wo2;
      z_den4[k] = - 4.0*a + 2.0*b - 2.0*b*wo2 + 4.0*a*wo2*wo2;
      z_den5[k] = a - b + 2.0*wo2*a + c -b*wo2 +a*wo2*wo2;
    }
    if (Nodd)
    {
      z_den1[kmax-1] = 1.0 + wc*sinh(vo) + wo2; //b + c + b*wo2;
      z_den2[kmax-1] = 2.0*(wo2-1.0);           //2.0*b*(wo2-1.0);
      z_den3[kmax-1] = 1.0 - wc*sinh(vo) + wo2; //b - c + b*wo2;
      z_den4[kmax-1] = 0.0;
      z_den5[kmax-1] = 0.0;
    }

    // Normalize denominator coeffs to ease filter implementation
    for (k=0;k<kmax;k++)
    {
      if ((z_den1[k] != 0.0) && (z_den1[k] != 1.0))
      {
        z_const = z_const/z_den1[k];
        z_den2[k] = z_den2[k]/z_den1[k];
        z_den3[k] = z_den3[k]/z_den1[k];
        z_den4[k] = z_den4[k]/z_den1[k];
        z_den5[k] = z_den5[k]/z_den1[k];
        z_den1[k] = 1.0;
      }
    }

    // initialize storage for delay elements, etc...
    w_del1 = new double[kmax]; // first delay element
    w_del2 = new double[kmax]; // second delay element
    w_del3 = new double[kmax]; // first delay element
    w_del4 = new double[kmax]; // second delay element
    w_delp = new double[kmax]; // present element (for iteration)
    for (k=0;k<kmax;k++)
    {
      w_del1[k] = 0.0;
      w_del2[k] = 0.0;
      w_del3[k] = 0.0;
      w_del4[k] = 0.0;
      w_delp[k] = 0.0;
    }
    // initialize intermediate storage for debug tracking...
    w = new double[kmax];
    x_int = new double[kmax];
    y_int = new double[kmax];
    for (k=0;k<kmax;k++)
    {
      w[k] = 0.0;
      x_int[k] = 0.0;
      y_int[k] = 0.0;
    }
  } // end of else block (bandpass filter initialization)

  // Debugging...
  //    if (!rep) { // put debugging info in std output only...
  //       cout << endl;
  //       cout << "****************************************************"
  //          << endl;
  //       cout << "Info for Zchebyi instance: " << name_ch;
  //       cout << endl;
  //       cout << endl;
  //       cout << "Parameters (and interpretation):" << endl;
  //       cout << "   fo    = " << fo << endl;
  //       cout << "   pbw   = " << pbw << endl;
  //       cout << "   pord  = " << pord << endl;
  //       cout << "   fls   = " << fls << endl;
  //       cout << "   flp   = " << flp << endl;
  //       cout << "   fhp   = " << fhp << endl;
  //       cout << "   fhs   = " << fhs << endl;
  //       cout << "   pbfdb = " << pbfdb << endl;
  //       cout << "   sbadb = " << sbadb << endl;
  //       cout << "   ildb  = " << ildb << endl;
  //       cout << endl;
  //
  //       cout << "Derived filter variables:" << endl;
  //       cout << "   Ngiven  = " << Ngiven << endl;
  //       cout << "   lowpass = " << lowpass << endl;
  //       cout << "   epsilon = " << epsilon << endl;
  //       cout << "   lambda  = " << lambda << endl;
  //       cout << "   tee     = " << tee << endl;
  //       cout << "   wls     = " << wls << endl;
  //       cout << "   wlp     = " << wlp << endl;
  //       cout << "   whp     = " << whp << endl;
  //       cout << "   whs     = " << whs << endl;
  //       cout << "   wp      = " << wp << endl;
  //       cout << "   ws      = " << ws << endl;
  //       cout << "   wo2     = " << wo2 << endl;
  //       cout << "   wo      = " << wo  << endl;
  //       cout << "   N       = " << N   << endl;
  //       cout << "   wc      = " << wc  << endl;
  //       cout << "   Nodd    = " << Nodd << endl;
  //       cout << "   kmax    = " << kmax << endl;
  //       cout << endl;
  //
  //       cout << "Lowpass Prototype specifications:" << endl;
  //       for (k=0;k<kmax;k++) {
  //       cout << "   plp_den(" << k << ") = " << plp_den1[k] << " "
  //          << plp_den2[k] << " " << plp_den3[k] << endl;
  //       }
  //       cout << "   plp_num    = " << plp_num << endl;
  //       cout << endl;
  //
  //       if (lowpass) {
  //          cout << "Lowpass Z-filter coefficients:" << endl;
  //          cout << "   z_const  = " << z_const << endl;
  //          for (k=0;k<kmax;k++) {
  //          cout << "   z_num(" << k << ") = " << z_num1[k] << " "
  //             << z_num2[k] << " " << z_num3[k] << endl;
  //          }
  //          for (k=0;k<kmax;k++) {
  //          cout << "   z_den(" << k << ") = " << z_den1[k] << " "
  //             << z_den2[k] << " " << z_den3[k] << endl;
  //          }
  //          cout << endl;
  //       }
  //       else { // bandpass
  //          cout << "Bandpass Z-filter coefficients:" << endl;
  //          cout << "   z_const  = " << z_const << endl;
  //          for (k=0;k<kmax;k++) {
  //          cout << "   z_num(" << k << ") = " << z_num1[k] << " "
  //             << z_num2[k] << " " << z_num3[k] << endl;
  //          }
  //          for (k=0;k<kmax;k++) {
  //          cout << "   z_den(" << k << ") = " << z_den1[k] << " "
  //             << z_den2[k] << " " << z_den3[k] << " "
  //             << z_den4[k] << " " << z_den5[k] << endl;
  //          }
  //          cout << endl;
  //       }
  //       cout << "****************************************************"
  //         << endl;
  //    }
  //
  // Record above into freeda report if requested at netlist
  if (rep)
  {
   newLine();
    sepLine();
    newLine();
    sprintf(msg,"Info for Zchebyi instance: %s",
        name_ch );
    report(LOGFILE, msg);
    newLine();
    sprintf(msg,"Parameters (and interpretation):");
    report(LOGFILE, msg);
    sprintf(msg,"   fo    = %g", fo);
    report(LOGFILE, msg);
    sprintf(msg,"   pbw   = %g", pbw);
    report(LOGFILE, msg);
    sprintf(msg,"   pord  = %i", pord);
    report(LOGFILE, msg);
    sprintf(msg,"   fls   = %g", fls);
    report(LOGFILE, msg);
    sprintf(msg,"   flp   = %g", flp);
    report(LOGFILE, msg);
    sprintf(msg,"   fhp   = %g", fhp);
    report(LOGFILE, msg);
    sprintf(msg,"   fhs   = %g", fhs);
    report(LOGFILE, msg);
    sprintf(msg,"   pbfdb = %g", pbfdb);
    report(LOGFILE, msg);
    sprintf(msg,"   sbadb = %g", sbadb);
    report(LOGFILE, msg);
    sprintf(msg,"   ildb  = %g", ildb);
    report(LOGFILE, msg);
    newLine();

    sprintf(msg,"Derived filter variables:");
    report(LOGFILE, msg);
    sprintf(msg,"   Ngiven  = %u", Ngiven);
    report(LOGFILE, msg);
    sprintf(msg,"   lowpass = %u", lowpass);
    report(LOGFILE, msg);
    sprintf(msg,"   epsilon = %g", epsilon);
    report(LOGFILE, msg);
    sprintf(msg,"   lambda  = %g", lambda);
    report(LOGFILE, msg);
    sprintf(msg,"   tee     = %g", tee);
    report(LOGFILE, msg);
    sprintf(msg,"   wls     = %g", wls);
    report(LOGFILE, msg);
    sprintf(msg,"   wlp     = %g", wlp);
    report(LOGFILE, msg);
    sprintf(msg,"   whp     = %g", whp);
    report(LOGFILE, msg);
    sprintf(msg,"   whs     = %g", whs);
    report(LOGFILE, msg);
    sprintf(msg,"   wp      = %g", wp);
    report(LOGFILE, msg);
    sprintf(msg,"   ws      = %g", ws);
    report(LOGFILE, msg);
    sprintf(msg,"   wo2     = %g", wo2);
    report(LOGFILE, msg);
    sprintf(msg,"   wo      = %g", wo);
    report(LOGFILE, msg);
    sprintf(msg,"   N       = %i", N);
    report(LOGFILE, msg);
    sprintf(msg,"   wc      = %g", wc);
    report(LOGFILE, msg);
    sprintf(msg,"   Nodd    = %u", Nodd);
    report(LOGFILE, msg);
    sprintf(msg,"   kmax    = %i", kmax);
    report(LOGFILE, msg);
    newLine();

    sprintf(msg,"Lowpass Prototype specifications:");
    report(LOGFILE, msg);
    for (k=0;k<kmax;k++)
    {
      sprintf(msg,"   plp_den(%i) = %g %g %g", k, plp_den1[k],
          plp_den2[k], plp_den3[k]);
      report(LOGFILE, msg);
    }
    sprintf(msg,"   plp_num    = %g", plp_num);
    report(LOGFILE, msg);
    newLine();

    if (lowpass)
    {
      sprintf(msg,"Lowpass Z-filter coefficients:");
      report(LOGFILE, msg);
      sprintf(msg,"   z_const  = %g", z_const);
      report(LOGFILE, msg);
      for (k=0;k<kmax;k++)
      {
        sprintf(msg,"   z_num(%i) = %g %g %g", k, z_num1[k],
            z_num2[k], z_num3[k]);
        report(LOGFILE, msg);
      }
      for (k=0;k<kmax;k++)
      {
        sprintf(msg,"   z_den(%i) = %g %g %g", k, z_den1[k],
            z_den2[k], z_den3[k]);
        report(LOGFILE, msg);
      }
      newLine();
    }
    else
    { // bandpass
      sprintf(msg,"Bandpass Z-filter coefficients:");
      report(LOGFILE, msg);
      sprintf(msg,"   z_const  = %g", z_const);
      report(LOGFILE, msg);
      for (k=0;k<kmax;k++)
      {
        sprintf(msg,"   z_num(%i) = %g %g %g", k, z_num1[k],
            z_num2[k], z_num3[k]);
        report(LOGFILE, msg);
      }
      for (k=0;k<kmax;k++)
      {
        sprintf(msg,"   z_den(%i) =  %g %g %g %g %g", k, z_den1[k],
            z_den2[k], z_den3[k], z_den4[k], z_den5[k]);
        report(LOGFILE, msg);
      }
      newLine();
    }
    sepLine();
    newLine();
  } // record to report()
}
