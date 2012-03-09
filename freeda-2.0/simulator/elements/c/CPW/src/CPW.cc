#include "../../../../network/Circuit.h"
#include "CPW.h"

// Static members
const unsigned CPW::n_par = 9;

// Element information
ItemInfo CPW::einfo =
{
  "cpw",
  "Coplanar transmission line",
  "Carlos E. Christoffersen",
  DEFAULT_ADDRESS"category:transmission line",
  "2000_07_20"
};

// Parameter information
ParmInfo CPW::pinfo[] =
{
  {"w", "Strip to ground space (m)", TR_DOUBLE, true},
  {"s", "Strip width (m)", TR_DOUBLE, true},
  {"length", "Physical length (m)", TR_DOUBLE, true},
  {"er", "Substrate dielectric constant", TR_DOUBLE, false},
  {"t", "Substrate thickness (m)", TR_DOUBLE, true},
  {"tand", "Loss tangent", TR_DOUBLE, false},
  {"sigma", "Metal conductivity (S/m)", TR_DOUBLE, false},
  {"nsect", "Enable discrete approximation with n sections", TR_INT, false},
  {"fopt", "Optimum frequency for discrete approximation", TR_DOUBLE, false}
};


CPW::CPW(const string& iname) : Element(&einfo, pinfo, n_par, iname)
{
  //
  paramvalue[0] = &w;
  paramvalue[1] = &s;
  paramvalue[2] = &length;
  paramvalue[3] = &(er = one);
  paramvalue[4] = &t;
  paramvalue[5] = &(tand = zero);
  paramvalue[6] = &(sigma = 5.8e7);
  paramvalue[7] = &(nsect = 0);
  paramvalue[8] = &(fopt = zero);

  // Set the number of terminals
  setNumTerms(4);

  setFlags(LINEAR | MULTI_REF | TR_FREQ_DOMAIN);

  tline = 0;
}

void CPW::init() throw(string&)
{
  Circuit* cir = getCircuit();
  unsigned newelem_id =
	cir->addElement("tlinp4", getInstanceName() + ":tline", true);
  tline = (Tlinp4*) cir->getElement(newelem_id);
  cir->connect(newelem_id, getTerminal(0)->getId());
  cir->connect(newelem_id, getTerminal(1)->getId());
  cir->connect(newelem_id, getTerminal(2)->getId());
  cir->connect(newelem_id, getTerminal(3)->getId());

  // Calculate the transmission line parameters
  // effective epsilon
  double ee = (er + one) / 2.;

  double k = s / (s + w * 2.);

  // The following K approximation is valid for k <= .4
  //
  // (Later add the formulas for other cases)
  if (k > .4)
    throw(getInstanceName() + ": k > .4 (add the formulas from Collin book)");
  double K = one + k*k/4. + k*k*k*k*9./64;
  K *= pi / 2.;

  assert(k <= .7);
  double kp = sqrt(one - k*k);
  double sqkp = sqrt(kp);
  // K/K' = KoKp
  double KoKp = pi / log((one + sqkp) * 2. / (one - sqkp));
  // characteristic impedance
  double z0 = sqrt(mu0 / epsilon0);
  double zc = z0 / 4. / sqrt(ee) / KoKp;

  // Now the attenuation
  // Calculate at fref = 1 GHz    att in (dB/m)
  const double fref = 1e9;
  double Rm = pi * sqrt(4e-7 / sigma) * sqrt(fref);
  // This is for Cu
  //  double Rm = 2.608e-7*sqrt(fref);
  double R1 = Rm / (4.*s*(one-k*k)*K*K);
  double R2 = R1;
  R1 *= (pi + log(4.*pi*s/t) - k * log((one+k)/(one-k)));
  R2 *= k * (pi + log(4.*pi*(s+2.*w)/t) - log((one+k)/(one-k)) / k);
  double alphac = (R1 + R2) / zc / 2.;

  // Set up transmission line parameters
  tline->setParam("k", &ee, TR_DOUBLE);
  tline->setParam("alpha", &alphac, TR_DOUBLE);
  tline->setParam("z0mag", &zc, TR_DOUBLE);
  tline->setParam("length", &length, TR_DOUBLE);
  tline->setParam("fscale", &fref, TR_DOUBLE);
  tline->setParam("tand", &tand, TR_DOUBLE);
  tline->setParam("nsect", &nsect, TR_INT);
  tline->setParam("fopt", &fopt, TR_DOUBLE);

  // init tline
  tline->init();
}

void CPW::getLocalRefIdx(UnsignedVector& local_ref_vec,
TerminalVector& term_list)
{
  tline->getLocalRefIdx(local_ref_vec, term_list);
}

void CPW::fillMNAM(FreqMNAM* mnam)
{
  return;
}

void CPW::fillMNAM(TimeMNAM* mnam)
{
  return;
}

