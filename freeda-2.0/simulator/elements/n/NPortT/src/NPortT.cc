#include "../../../../network/ElementManager.h"
#include "../../../../analysis/FreqMNAM.h"
#include "../../../../analysis/TimeMNAM.h"
#include "NPortT.h"

#include <cstdio>
#include <cstdlib>

extern "C"
{
	#include "../../../../compat/in_interp.h"
}
#include "../../../lex.NPortT_.c"

// Static members
const unsigned NPortT::n_par = 2;

// Element information
ItemInfo NPortT::einfo =
{
  "nportt",
  "Multi-port element defined by its y parameter transfer functions",
  "Carlos E. Christoffersen",
  DEFAULT_ADDRESS"category:multiport,electrothermal",
  "2000_07_20"
};

// Parameter information
ParmInfo NPortT::pinfo[] =
{
  {"filename",
	"File containing the port-based parameter matrix.", TR_STRING, true},
  {"mf", "Number to multiply the normalizing factor (time domain only)",
	TR_DOUBLE, false}
};


NPortT::NPortT(const string& iname) : Element(&einfo, pinfo, n_par, iname)
{
  // Standard construction code.
  paramvalue[0] = &filename;
  paramvalue[1] = &(mf = one);

  setFlags(LINEAR | MULTI_REF | TR_FREQ_DOMAIN);

  // Set initial terminal and reference count from file.
  ports = 0;
  current_group = 0;

  my_start_row = 0;
  extra_rcs = 0;
  denvec = numvec = 0;
}


NPortT::~NPortT()
{
  if (denvec) {
    delete [] denvec;
    delete [] numvec;
  }
}

void NPortT::init() throw(string&)
{
  FILE* fp;
  // Open input file
  if( (fp = fopen(filename.c_str(), "r")) == NULL)
    throw(getInstanceName() + ": Cannot open input file \"" +
	filename +"\".");
  // Point scanner pointer here.
  scan_NPortT.thisNPortT = this;
  NPortT_in = fp;
  // Set initial value for counters
  scan_NPortT.ncoeff = 0;
  scan_NPortT.coeffcount = 0;
  scan_NPortT.trcount = 0;
  scan_NPortT.linecount = 1;
  // Parse first part (ports) of input file
  if(NPortT_lex())
    throw(getInstanceName() + ": in file \"" + filename + "\": " +
	+ scan_NPortT.serror);
  // Check that the number of terminals in netlist is consistent
  // with number of terminals in file
  if (getTermCount() != ports + group_vec.size())
    throw(string(getInstanceName() + ": Incorrect number of terminals."));
  // If control reaches here, then the number of terminals in
  // netlist is correct.
  setNumTerms(getTermCount());
  // Set the number of transfer functions
  ntransfer = ports * ports;
  // Allocate space for transfer functions
  numvec = new DenseDoubleVector[ntransfer];
	for (int i = 0; i < ntransfer; i++)
		numvec[i] = DenseDoubleVector(0);
  denvec = new DenseDoubleVector[ntransfer];
	for (int i = 0; i < ntransfer; i++)
		denvec[i] = DenseDoubleVector(0);
  // resize workspace matrix
  ws_matrix.reshape(getNumTerms(), getNumTerms());
  // resize auxiliary vectors
  neqns_row.resize(ports);
  neqns_col.resize(ports);
  // Continue parsing the rest of the file (y parameters).
  if (NPortT_lex())
    throw(getInstanceName() + ": in file \"" + filename + "\": " +
	+ scan_NPortT.serror);
  // Close input file
  fclose(fp);
  assert(denvec[ntransfer-1].length() > 0);
}

void NPortT::addPort(const unsigned& group_number)
{
  if (!ports || current_group != group_number)
	{
    // This is the first terminal or
    // there is a new group.
    // Add a new group with one count (the current terminal).
    group_vec.push_back(1);
    current_group = group_number;
  }
  else
    // Increment the terminal counter for this group.
	(group_vec.back())++;

  // Increment the number of ports given in the file.
  ports++;
}

void NPortT::setNumSize(unsigned tr_idx, unsigned size)
{
  assert(tr_idx < ntransfer);
  numvec[tr_idx].resize(size);
}

void NPortT::setDenSize(unsigned tr_idx, unsigned size)
{
  assert(tr_idx < ntransfer);
  denvec[tr_idx].resize(size);
}

double& NPortT::setNumCoeff(unsigned tr_idx, unsigned idx)
{
  assert(tr_idx < ntransfer);
  return numvec[tr_idx][idx];
}

double& NPortT::setDenCoeff(unsigned tr_idx, unsigned idx)
{
  assert(tr_idx < ntransfer);
  return denvec[tr_idx][idx];
}

void NPortT::getLocalRefIdx(UnsignedVector& local_ref_vec,
TerminalVector& term_list)
{
  // Make sure the vectors are empty
  term_list.erase(term_list.begin(), term_list.end());
  local_ref_vec.erase(local_ref_vec.begin(), local_ref_vec.end());

  unsigned termcount = 0;
  // Iterate for each group
  for (unsigned i = 0; i < group_vec.size(); i++)
	{
    // Iterate for each port in group
    for (unsigned j = 0; j < group_vec[i]; j++)
		{
      term_list.push_back(getTerminal(termcount));
      // Go to next terminal (in netlist connection).
      termcount++;
    }
    // Add now the local reference terminal
    term_list.push_back(getTerminal(ports + i));
    local_ref_vec.push_back(termcount + i); // Local reference index
  }
}

unsigned NPortT::getExtraRC(const unsigned& eqn_number, const MNAMType& type)
{
  if (type == TIME_DOMAIN)
	{
    // Find the maximun degree of the numerator in each row -> current eqns.
    ncurr_eqns = 0;
    for (unsigned i = 0; i < ports; i++)
		{
      neqns_row[i] = 0;
      for (unsigned j = 0; j < ports; j++)
			{
				// get the degree of the poly
				unsigned ntmp = denvec[i + j * ports].length() - 1;
				if (ntmp > unsigned(neqns_row[i]))
					neqns_row[i] = ntmp;
      }
      ncurr_eqns += neqns_row[i];
    }
    // Find the maximun degree of the denominator
    // in each column -> voltage eqns.
    nvolt_eqns = 0;
    for (unsigned i = 0; i < ports; i++)
		{
      neqns_col[i] = 0;
      for (unsigned j = 0; j < ports; j++)
			{
				// get the degree of the poly
				unsigned ntmp = numvec[i * ports + j].length() - 1;
				if (ntmp > unsigned(neqns_col[i]))
					neqns_col[i] = ntmp;
      }
      nvolt_eqns += neqns_col[i];
    }
    // Keep the equation number assigned to this element
    my_start_row = eqn_number;
    // Add extra RCs
    extra_rcs = ncurr_eqns + nvolt_eqns + ports;
  }
  else
	{
    my_start_row = 0;
    extra_rcs = 0;
  }
  return extra_rcs;
}

void NPortT::getExtraRC(unsigned& first_eqn, unsigned& n_rows) const
{
  first_eqn = my_start_row;
  n_rows = extra_rcs;
}

void NPortT::fillMNAM(FreqMNAM* mnam)
{
  double_complex jw = mnam->getFreq() * double_complex(0., 2.) / fmax;
  // Do not attempt to use the model out of range
  if (mnam->getFreq() > fmax)
    jw = fmax * double_complex(0., 2.);

  for (unsigned i = 0; i < ports; i++)
    for (unsigned j = 0; j < ports; j++)
	{
		double_complex num = numvec[i + j * ports][0];
		unsigned kmax = numvec[i + j * ports].length();
		for (unsigned k = 1; k < kmax; k++)
		{
			num *= jw;
			num += numvec[i + j * ports][k];
		}
		double_complex den = denvec[i + j * ports][0];
		kmax = denvec[i + j * ports].length();
		for (unsigned k = 1; k < kmax; k++)
		{
			den *= jw;
			den += denvec[i + j * ports][k];
		}
    ws_matrix(i,j) = (num/den);
	}
	// Now generate the missing rows and columns.
	unsigned lbase = 0;
	for (unsigned k = 0; k < group_vec.size(); k++)
	{
		// Iterate for each port in group
		for (unsigned i = 0; i < ports; i++)
		{
			for (unsigned l = 0; l < group_vec[k]; l++)
			{
        ws_matrix(ports+k,i) -= ws_matrix(l+lbase,i);
        ws_matrix(i,ports+k) -= ws_matrix(i,l+lbase);
			}
		}
		lbase += group_vec[k];
	}
	lbase = 0;
	for (unsigned k = 0; k < group_vec.size(); k++)
	{
		// Iterate for each port in group
		for (unsigned i = ports; i < getNumTerms(); i++)
			for (unsigned l = 0; l < group_vec[k]; l++)
        ws_matrix(i,ports+k) -= ws_matrix(i,l+lbase);
		lbase += group_vec[k];
	}
	// Finally, fill the MNAM
	for (unsigned i = 0; i < getNumTerms(); i++)
		for (unsigned j = 0; j < getNumTerms(); j++)
			mnam->setElement(getTerminal(i)->getRC(),
	getTerminal(j)->getRC(),
	ws_matrix(i,j));
}

void NPortT::fillMNAM(TimeMNAM* mnam)
{
  // Remember that all this is valid only if p0 is normalized to 1
  // where H(s) = Q(s) / P(s)
  // Q(s) = q0 + q1 s + q2 s^2 + ... + qn s^n
  // P(s) = p0 + p1 s + p2 s^2 + ... + pv s^v

  // ********************************************************************
  // Fill current equations
  // ********************************************************************
  // Base for this row's port current derivatives
  unsigned ii_counter = my_start_row + ports;
  // i: row index
  unsigned i = 0;
  for (unsigned m = 0; m < group_vec.size(); m++)
	{
    for (unsigned n = 0; n < group_vec[m]; n++, i++)
		{
      // this current equation number
      unsigned ceqnumb = my_start_row + i;
      mnam->setMElement(getTerminal(i)->getRC(), ceqnumb, one);
      mnam->setMElement(getTerminal(ports + m)->getRC(), ceqnumb, -one);
      mnam->setMElement(ceqnumb, ceqnumb, -one);
      // Base for this column's port voltage derivatives
      unsigned vj_counter = my_start_row + ports + ncurr_eqns;
      // j: column index
      unsigned j = 0;
      for (unsigned k = 0; k < group_vec.size(); k++)
			{
				for (unsigned l = 0; l < group_vec[k]; l++, j++)
				{
					DenseDoubleVector& num = numvec[i + j * ports];
					const unsigned& numsize = num.length();
					// fill current equation
					if (num[numsize - 1])
					{
						mnam->setMElement(ceqnumb, getTerminal(j)->getRC(),
			      num[numsize - 1]);
						mnam->setMElement(ceqnumb, getTerminal(ports + k)->getRC(),
			      -num[numsize - 1]);
					}
					for (unsigned vcount = 0; vcount < numsize - 1; vcount++)
						if (num[numsize - vcount - 2])
							mnam->setMElement(ceqnumb, vj_counter + vcount,
					num[numsize - vcount - 2] / mf);
					DenseDoubleVector& den = denvec[i + j * ports];
					const unsigned& densize = den.length();
					for (unsigned icount = 0; icount < densize - 1; icount++)
						if (den[densize - icount - 2])
							mnam->setMElement(ceqnumb, ii_counter + icount,
					-den[densize - icount - 2] / mf);
					// Increment base for this column's port voltage derivatives
					vj_counter += neqns_col[j];
				}
      }
      // Increment base for this row's port current derivatives
      ii_counter += neqns_row[i];
    }
  }

  // ********************************************************************
  // Fill current and voltage derivative equations
  // ********************************************************************
  double factor = one / fmax / pi;
  // Base for this row's port current derivatives
  ii_counter = my_start_row + ports;
  for (unsigned i = 0; i < ports; i++)
	{
    unsigned ceqnumb = my_start_row + i;
    mnam->setMpElement(ii_counter, ceqnumb, - mf * factor);
    mnam->setMElement(ii_counter, ii_counter, one);
    for (unsigned icount = 1; icount < unsigned(neqns_row[i]); icount++)
		{
      mnam->setMpElement(ii_counter + icount, ii_counter + icount - 1,
			- mf * factor);
      mnam->setMElement(ii_counter + icount, ii_counter + icount, one);
    }
    // Increment base for this row's port current derivatives
    ii_counter += neqns_row[i];
  }

  // Base for this row's port voltage derivatives
  ii_counter = my_start_row + ports + ncurr_eqns;
  i = 0;
  for (unsigned m = 0; m < group_vec.size(); m++)
	{
    for (unsigned n = 0; n < group_vec[m]; n++, i++)
		{
      mnam->setMpElement(ii_counter, getTerminal(i)->getRC(), - mf * factor);
      mnam->setMpElement(ii_counter, getTerminal(ports + m)->getRC(),
			mf * factor);
      mnam->setMElement(ii_counter, ii_counter, one);
      for (unsigned icount = 1; icount < unsigned(neqns_col[i]); icount++)
			{
				mnam->setMpElement(ii_counter + icount, ii_counter + icount - 1,
				- mf * factor);
				mnam->setMElement(ii_counter + icount, ii_counter + icount, one);
      }
      // Increment base for this row's port voltage derivatives
      ii_counter += neqns_col[i];
    }
  }
}

