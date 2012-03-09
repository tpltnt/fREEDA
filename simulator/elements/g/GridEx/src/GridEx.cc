#include "../../../../network/ElementManager.h"
#include "../../../../analysis/FreqMNAM.h"
#include "../../../../analysis/TimeMNAM.h"
#include "../../../../analysis/TimeDomainSV.h"
#include "GridEx.h"

#include <cstdlib>

extern "C" {
#include "../../../../compat/in_interp.h"
	   }

#define MAX_SPICE_LINE_LEN 160

// Static members
const unsigned GridEx::n_par = 4;

// Element information
ItemInfo GridEx::einfo =
{
  "gridex",
  "Electromagnetic grid excitation",
  "Mark Summers, Carlos E. Christoffersen",
  DEFAULT_ADDRESS"category:multiport",
  "2000_09_10"
};

// Parameter information
ParmInfo GridEx::pinfo[] =
{
  {"ifilename",
   "File containing normalized currents.", TR_STRING, true},
  {"efilename",
   "File containing field coefficients for currents.", TR_STRING, true},
  {"freq",
   "Frequency of the excitation (Hz).", TR_DOUBLE, false},
  {"delay", "Time before which the output current is set to 0 (s)",
   TR_DOUBLE, false}
};


GridEx::GridEx(const string& iname) : Element(&einfo, pinfo, n_par, iname)
{
  // Standard construction code.
  paramvalue[0] = &ifilename;
  paramvalue[1] = &efilename;
  paramvalue[2] = &(my_freq = zero);
  paramvalue[3] = &(delay = zero);


  setFlags(LINEAR | MULTI_REF | TR_TIME_DOMAIN | SOURCE);

  // For checking purposes
  ivec_size = 0;
  ifreq_vec = NULL;
  i_vec = NULL;

  evec_size = 0;
  efreq_vec = NULL;
  e_vec = NULL;

  last_i_index = last_e_index = 0;
}


GridEx::~GridEx()
{
  // Free allocated vectors
  if (ivec_size) {
    delete [] ifreq_vec;
    delete [] i_vec;
  }
  if (evec_size) {
    delete [] efreq_vec;
    delete [] e_vec;
  }
}


void GridEx::check() throw(string&)
{
  // We need to check parameters first, before attepmting to open
  // the data file.
  try {
    checkParams();
  }
  catch(string& msg) {
    throw(getInstanceName() + ": " + msg);
  }
}


void GridEx::init() throw(string&)
{
  FILE* fp;
  int f_count;

  // Open current input file
  if( (fp = fopen(ifilename.c_str(), "r")) == NULL)
    throw(getInstanceName() + ": Cannot open input file \"" +
	  ifilename +"\".");

  if (getTermCount() % 2)
    throw(getInstanceName() +
	  ": incorrect number of terminals (must be even)");

  // Check that the number of terminals is correct
  f_count = iac_find_count(fp);
  if (f_count % (getTermCount() + 1))
    throw(getInstanceName() + ": wrong number of values in data file \"" +
	  ifilename +"\".");

  ivec_size = f_count / (getTermCount() + 1);
  // Now we know the number of terminals is correct
  setNumTerms(getTermCount());
  nsources = getNumTerms()>>1;

  // Allocate appropiate amount of memory
  ifreq_vec = new double[ivec_size];
  i_vec = new DenseComplexVector[ivec_size];
	for (int i = 0; i < ivec_size; i++)
		i_vec[i] = DenseComplexVector(nsources);

  /*  Goto begining of file  */
  rewind(fp);

  /* store data from file */
  iac_store(i_vec, ifreq_vec, fp, ivec_size, nsources);

  /* close file */
  fclose(fp);


  /*  Now get E-field Normalizing Constants */
  // Open field input file
  if( (fp = fopen(efilename.c_str(), "r")) == NULL)
    throw(getInstanceName() + ": Cannot open input file \"" +
	  efilename +"\".");

  // Check number of values in files
  f_count = iac_find_count(fp);
  if(f_count % 3)
    throw(getInstanceName() + ": wrong number of values in data file \"" +
	  efilename +"\".");
  // Set size of field frequency vector
  evec_size = f_count / 3;

  efreq_vec = new double[evec_size];
  e_vec = new DenseComplexVector[evec_size];
	for (int i = 0; i < evec_size; i++)
		e_vec[i] = DenseComplexVector(1);

  /*  Goto begining of file  */
  rewind(fp);

  /* store data from file */
  iac_store(e_vec, efreq_vec, fp, evec_size, 1);

  fclose(fp);	 /* close file */

  // Set number of states
  setNumberOfStates(nsources);

  // Now init variables for transient
  // Get field
  double_complex e_val = get_field(my_freq);

  tr_ivec.resize(nsources);
  for (unsigned i=0; i<nsources; i++)
    tr_ivec[i] = get_current(my_freq, i) * e_val;
}

void GridEx::getLocalRefIdx(UnsignedVector& local_ref_vec,
			     TerminalVector& term_list)
{
  // Make sure the vectors are empty
  term_list.erase(term_list.begin(), term_list.end());
  local_ref_vec.erase(local_ref_vec.begin(), local_ref_vec.end());

  // The terminal order is alternating terminal and reference
  for (unsigned i = 0; i < getNumTerms(); i += 2) {
    term_list.push_back(getTerminal(i));
    term_list.push_back(getTerminal(i+1)); // Local reference
    local_ref_vec.push_back(i+1); // Local reference index
  }
}


void GridEx::fillMNAM(FreqMNAM* mnam)
{
  const double& freq(mnam->getFreq());
  // Fill source vector only at the desired frequency
  if (!isSet(&my_freq) || (freq && abs(freq - my_freq)/freq < 1e-14 )) {

    double_complex i_val, e_val;

    // Get the field value for this frequency.
    e_val = get_field(freq);

    // Fill the source vector for each port
    for (unsigned j=0; j < nsources; j++) {
      i_val = get_current(freq, j) * e_val;
      mnam->addToSource(getTerminal(2*j)->getRC(),
			getTerminal(2*j+1)->getRC(),
			i_val);
    }
  }
}

void GridEx::fillMNAM(TimeMNAM* mnam)
{
  // Nothing to be done.
  return;
}

void GridEx::fillSourceV(TimeMNAM* mnam)
{
  const double& ctime = mnam->getTime();
  // Calculate each current value
  if (ctime >= delay)
    for (unsigned j=0; j < nsources; j++) {
      double i_val = abs(tr_ivec[j]) *
	cos(twopi * my_freq * ctime + arg(tr_ivec[j]));
      mnam->addToSource(getTerminal(2*j)->getRC(),
			getTerminal(2*j+1)->getRC(),
			i_val);
    }
}

void GridEx::svTran(TimeDomainSV* tdsv)
{
  const double& ctime = tdsv->getCurrentTime();
  // Calculate each current value
  if (ctime < delay || tdsv->DC())
    for (unsigned j=0; j < nsources; j++) {
      tdsv->i(j) = zero;
      tdsv->u(j) = tdsv->getX(j);
    }
  else
    for (unsigned j=0; j < nsources; j++) {
      tdsv->i(j) =  - abs(tr_ivec[j]) *
	cos(twopi * my_freq * ctime + arg(tr_ivec[j]));
      tdsv->u(j) = tdsv->getX(j);
    }
}

void GridEx::deriv_svTran(TimeDomainSV* tdsv)
{
  tdsv->cleanJac();
  DoubleDenseMatrix& Ju = tdsv->getJu();
  for (unsigned j=0; j < nsources; j++) {
    Ju(j,j) = one;
  }
}


int GridEx::iac_find_count(FILE *fp)
{
  char c;
  char line[MAX_SPICE_LINE_LEN];
  int tmp_count = 0;


  while ((c = fgetc(fp)) != EOF) {
    if(!isspace(c)) {
      switch(c) {
	/* Indicates comment.  Ignore */
      case '!':
	fgets(line, MAX_SPICE_LINE_LEN, fp);
	break;

	/* Look for digit, decimal point or negative sign to
	 * signal start of number */

      case '0': case '1': case '2': case '3': case '4':  case '5':
      case '6': case '7': case '8': case '9': case '-': case '.':

	/* Look for white space that signals end of number */
	while(!isspace(c)) {
	  if(feof(fp)) {
		break;
	      }
	  c = fgetc(fp);
	}
	/*  Increment count for each number found  */
	tmp_count++;
	break;
      default:
	break;
      } /* Switch  */
    } /*  If  */
  } /* while */

  return(tmp_count);
}



int GridEx::iac_store(DenseComplexVector* c_vec_P,
		      double *f_vec, FILE *fp, int f_count, int t_count)
{
  char c;
  char line[MAX_SPICE_LINE_LEN];
  int f,j,i;
  double re, im;

  while ((c = fgetc(fp)) != EOF) {
    if(!isspace(c)) {
      switch(c) {
      case '!':
	fgets(line, MAX_SPICE_LINE_LEN, fp);
	break;

      default:

	ungetc(c,fp); /* Push the last character read back on stream */

	/*  Iterate thru the file and store the data using fscanf()  */

	for(f = 0; f< f_count; f++) {
	  i= fscanf(fp, "%lg", &(f_vec[f]));  /* Read in frequencies  */
	  while(i==0) {
	    if(fgetc(fp) == '!') {
	      fgets(line, MAX_SPICE_LINE_LEN, fp);
	    }
	    i= fscanf(fp, "%lg", &(f_vec[f]));
	  }
	  /* Now find the current information */
	  for(j=0; j<t_count; j++) {
	    i = fscanf(fp, "%lg", &(re));
	    /*  Remove any stray characters */
	    while(i==0) {
	      if(fgetc(fp) == '!') {
		fgets(line, MAX_SPICE_LINE_LEN, fp);
	      }
	      i = fscanf(fp, "%lg", &(re));
	    }
	    i = fscanf(fp, "%lg", &(im));
	    /*  Remove any stray characters */
	    while(i==0) {
	      if(fgetc(fp) == '!') {
		fgets(line, MAX_SPICE_LINE_LEN, fp);
	      }
	      i = fscanf(fp, "%lg", &(im));
	    }
	    assert(f < f_count);
	    assert(j < t_count);
	    (c_vec_P[f])[j] = double_complex(re,im);

	  } /*j*/
	} /* m */
      } /* Switch  */
    } /*  If  */
  } /* while */
  return 0;
}


double_complex GridEx::get_current(double f, unsigned port)
{
  double_complex val;
  assert(port < nsources);
  doubleHunt(ifreq_vec, ivec_size, f, &last_i_index);

  if(last_i_index <= 0) {
    val = i_vec[0][port];
  }
  else if(last_i_index >= int(ivec_size)) {
    val = i_vec[ivec_size-1][port];
  }
  else {
    // Interpolate
    val = i_vec[last_i_index][port] - i_vec[last_i_index-1][port];
    val *= ((f - ifreq_vec[last_i_index-1]) /
	    (ifreq_vec[last_i_index] - ifreq_vec[last_i_index-1]));
    val += i_vec[last_i_index-1][port];
  }
  return val;
}


double_complex GridEx::get_field(double f)
{
  double_complex val;
  doubleHunt(efreq_vec, evec_size, f, &last_e_index);

  if(last_e_index <= 0) {
    val = e_vec[0][0];
  }
  else if(last_e_index >= int(evec_size)) {
    val = e_vec[evec_size-1][0];
  }
  else {
    // Interpolate
    val = e_vec[last_e_index][0] - e_vec[last_e_index-1][0];
    val *= ((f - efreq_vec[last_e_index-1]) /
	    (efreq_vec[last_e_index] - efreq_vec[last_e_index-1]));
    val += e_vec[last_e_index-1][0];
  }
  return val;
}

