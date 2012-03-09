#include "../../../../network/ElementManager.h"
#include "../../../../analysis/FreqMNAM.h"
#include "Tlinp4.h"
#include <cstdio>

extern "C"
{
#include "../../../../inout/report.h"
#include "../../../../compat/standard.h"
}

// Tlinp
//
// This element models a physical transmission line.
// 
// Netlist form
// tlinp4:1 n1 n2 n3 n4 <parameter list>
// n2 and n4 are reference terminals
// 
// Two models are supported dependent on the secting of nsect.
//
// When nsect=0 (not set)
// In one model z0 and gamma are used in frequency domain analysis. This is
// selected with nsect is 0.
//
// When nsect>0
// The second model form expands the transmission line in a number of
// RLCG subsections.
// Terminal n4 should be the same as n2. The element then has one local
// reference terminal.
// This model can be used in the time or frequency domains.
//
// The model does not do the parameter checking it should.

// Static members

// Set the number of parameters. This number should be equal to the
// number of parametrers in ParmInfo below. The compiler is unable
// to automatically determine the length of the array as there is a
// variable length element (the string) embedded in the array.
const unsigned Tlinp4::n_par = 8;

// Set the element information. First entry should be in lower case and
// defines the name of the element model. The rets is used in self
// documentation.
ItemInfo Tlinp4::einfo =
{
  "tlinp4",
  "4 terminal physical transmission line",
  "Carlos E. Christoffersen, Mete Ozkar, Michael Steer",
  "category:transmission line",
  "2008_04_05"
};

// Parameter information
ParmInfo Tlinp4::pinfo[] =
{
  {"k", "Effective dielectric constant", TR_DOUBLE, false},
  {"alpha", "Attenuation (dB/m)", TR_DOUBLE, false},
  {"z0mag", "Magnitude of characteristic impedance (ohms)", TR_DOUBLE, true},
  {"fscale", "Scaling frequency for attenuation (Hz)", TR_DOUBLE, false},
  {"tand", "Loss tangent", TR_DOUBLE, false},
  {"length", "Line length (m)", TR_DOUBLE, true},
  {"nsect", "Enable discrete approximation with n sections", TR_INT, false},
  {"fopt", "Optimum frequency for discrete approximation", TR_DOUBLE, false}
};

// The constructor. Called when creating a new element instance.
Tlinp4::Tlinp4(const string& iname) : Element(&einfo, pinfo, n_par, iname)
{
  // Parameters are stored in the paraminfo vector but he we set default
  // values and the working name of the parameter.
  paramvalue[0] = &(k = one);
  paramvalue[1] = &(alpha = .1);
  paramvalue[2] = &z0mag;
  paramvalue[3] = &(fscale = zero);
  paramvalue[4] = &(tand = zero);
  paramvalue[5] = &length;
  paramvalue[6] = &(nsect = 0);
  paramvalue[7] = &(fopt = zero);

  // Set the number of terminals
  setNumTerms(4);

  // Set flags to tell the analysis routines what is special about this
  // element. (Here we are setting bits.)
  // LINEAR = this is a linear element
  // MULTIREF = this element has two or more reference terminals
  // TR_FREQ_DOMAIN = this element ed in transient or frequency domain analysis.
  // Both LINEAR and TR_FREQ_DOMAIN affect the type of analysis to be used.
  setFlags(LINEAR | MULTI_REF | TR_FREQ_DOMAIN);
}

// The initialization routine. Called prior to an analysis. Here the circuit
// can be changed and flags reset.
// The concept here is (if nsect is set so this model uses multiple sections)
// to build a circuit up of L's and Called's. The each L and C is 
// initialized individually as they will not be called otherwise.
void Tlinp4::init() throw(string&)
{
  // Calculate the capacitance per unit length
  c = sqrt(k)/(z0mag * c0);
  // Calculate the inductance per unit length
  l = (z0mag * z0mag) * c;
  // convert alpha from db/m to nepers/m
  alpha_nepers = alpha * 0.11512925;
  // If nsect (number of sections) is set then use a multi-section model for
  // the transmisison line.
  // The strategy is to expand the transmission line into a number of
  // RLGC sections. Attenuation circuit is built up with series L's and shunt C's
  // (which supported series resistance and shunt capacitance)
  if (nsect)
	{
    // Clear flags so this element is not called to fill the MNAM
    // Instead the individual L's and C's routines are used.
    setFlags(MULTI_REF); // So the element will not be used
                         // in linear and nonlinear analysis.

    // Need to check that the local reference terminals are the
    // same as this is required by the sectional model.
    if(getTerminal(1) != getTerminal(3))
    {
      char *cname;
      cname = (char *)(getInstanceName().c_str()); //get name of element
      char s[MAX_STRING_LEN+100];
      sprintf(s, "\nTerminal error for element %s . Parameter nsect > 0 and \
        sectional model requires element reference terminals to be the same.\n", cname);
      report(FATAL, s); // fREEDA will terminate with error message
    }

    // Use discrete approximation. Find the number of subsections and the
    // RLCG parameters of each.
    // First get the length of each subsection.
    double delta_x = length / double(nsect);
    double alphaf = (fscale == zero || fopt == zero) ? alpha_nepers :
      alpha_nepers * sqrt(fopt / fscale);
    // Now get the R, L, G and C of each subsection.
    double R = 2.0 * alphaf * z0mag * delta_x;
    double G = (fopt && tand) ? tand * twopi * fopt * c * delta_x : 1e-10;

    double L = l * delta_x;
    double C = c * delta_x;

    // Expand into nsect sections
    Circuit* cir = getCircuit();
    unsigned term_id1 = getTerminal(0)->getId();
    unsigned tref_id = getTerminal(1)->getId();
    for (int i=0; i<nsect; i++)
		{
      char loopit[10]; // so we can support 1,000,000,000 sections.
      // create a unique identifier
      sprintf(loopit, "%d", i);
      // Add serial branch (inductor with nopropagate flag)
      unsigned newelem_id = cir->addElement("inductor", 
        getInstanceName() + ":inductor:" + loopit, true);
      // Connect to previous terminal
      cir->connect(newelem_id, term_id1);
      // create new terminal (with nocheck) unless this is the last one
      unsigned term_id2;
      if (i != nsect-1)
        term_id2 = cir->addTerminal(getInstanceName() + ":" + loopit, true);
      else
        term_id2 = getTerminal(2)->getId();
      // Connect to term_id2
      cir->connect(newelem_id, term_id2);
      // Get inductor pointer.
      Element* elem = cir->getElement(newelem_id);
      // Set serial parameters
      elem->setParam("l", &L, TR_DOUBLE);
      elem->setParam("int_res", &R, TR_DOUBLE);
      // Initialize the inductor
      elem->init();

      // Add parallel branch (capacitor with nopropagate flag)
      newelem_id =
      cir->addElement("capacitor", 
        getInstanceName() + ":capacitor:" + loopit, true);
      // Connect to previous terminal and reference
      cir->connect(newelem_id, term_id2);
      cir->connect(newelem_id, tref_id);
      // Get capacitor pointer.
      elem = cir->getElement(newelem_id);
      // Set parallel parameters
      elem->setParam("c", &C, TR_DOUBLE);
      elem->setParam("int_g", &G, TR_DOUBLE);
      // Initialize the capacitor
      elem->init();

      // Prepare for next loop
      term_id1 = term_id2;
    }
  }
}

void Tlinp4::getLocalRefIdx(UnsignedVector& local_ref_vec,
			    TerminalVector& term_list)
{
  // Make sure the vectors are empty
  term_list.erase(term_list.begin(), term_list.end());
  local_ref_vec.erase(local_ref_vec.begin(), local_ref_vec.end());

  // Insert the four terminals in the vector listing terminals
  term_list.push_back(getTerminal(0));
  term_list.push_back(getTerminal(1)); // Local reference terminal
  term_list.push_back(getTerminal(2));
  term_list.push_back(getTerminal(3)); // Local reference terminal
  // Insert the two reference terminals in the vector listing reference
  // terminals of elements.
  local_ref_vec.push_back(1); // Local reference index
  local_ref_vec.push_back(3); // Local reference index
}


// Routine to fill the MNAM in frequency domain analysis.
void Tlinp4::fillMNAM(FreqMNAM* mnam)
{
  // If nsect is set, there is nothing to do
  // This is because the sectional model is being used and the MNAM will
  // be set by fillMNAM for the elements in the sectional model.
  if (nsect)
    return;

  // So we are not using the sectional model so use Z0 and gamma
  const double& freq(mnam->getFreq());
  // To avoid division by zero
  const double alphaMin = 1e-5;
  double omega = twopi * freq;

  double_complex Z0, gamma;
  if (freq == zero)
	{
    Z0 = z0mag;
    gamma = alpha_nepers;
  }
  else
	{
    double alphaf = (fscale == zero) ? alpha_nepers :
      alpha_nepers * sqrt(freq / fscale);
    double r = 2.0 * alphaf * z0mag;
    double g = tand * omega * c;

    double_complex z1(r , omega * l);
    double_complex y1(g , omega * c);
    // characterestic impedance
    // Zo = sqrt( (R + jwL)/(G + jwC))
    Z0 = sqrt(z1 / y1);
    // attenuation
    // gamma = sqrt( (R + jwL)*(G + jwC))
    gamma = sqrt(z1 * y1);
  }
  // Make sure that the real part of gamma is not too small
  if (gamma.real() < alphaMin)
    gamma = double_complex(alphaMin , gamma.imag());

  double_complex ctemp(gamma * length);

  double_complex y11(one);
  y11 /= Z0;
  y11 *= cosh(ctemp);
  y11 /= sinh(ctemp);

  double_complex y12(one);
  y12 /= Z0;
  y12 /= - sinh(ctemp);

  // Set admittances
  mnam->setAdmittance(getTerminal(0)->getRC(), getTerminal(1)->getRC(), y11);
  mnam->setAdmittance(getTerminal(2)->getRC(), getTerminal(3)->getRC(), y11);

  // Set voltage-controlled current sources
  mnam->setQuad(getTerminal(0)->getRC(), getTerminal(1)->getRC(),
		getTerminal(2)->getRC(), getTerminal(3)->getRC(),
		y12);
  mnam->setQuad(getTerminal(2)->getRC(), getTerminal(3)->getRC(),
		getTerminal(0)->getRC(), getTerminal(1)->getRC(),
		y12);
}

// Now we do not need to provide the time domain MNAM as we must have
// a sectional model for trasient analysis.
// If this routine is called and nsect = 0 then the sectional model has not
// been created and so there is no point continuing.
// If some types of transient analysis (such as a transoient analysis using
// convolution, the time domain MNAM is not required so sometimes it is
// valid to use nsect = 0. So it is only when this routine is called is
// nsect = 0 a problem.
void Tlinp4::fillMNAM(TimeMNAM* mnam)
{
  // The time-domain MNAM is set by the expanded circuit. However
  // if nsect =0 then there is not an expanded circuit
  if(nsect == 0)
  {
    char *cname;
    cname = (char *)(getInstanceName().c_str()); //get name of element instance
    char s[MAX_STRING_LEN+100];
    sprintf(s,
        "\nParameter error for element %s . Parametersrameter nsect = 0.\n\
        There are no sections.\n", cname);
    report(FATAL, s); // fREEDA will terminate with error message
  }
  assert(nsect);
}

