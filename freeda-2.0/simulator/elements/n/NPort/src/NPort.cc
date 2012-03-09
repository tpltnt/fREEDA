#include "../../../../network/ElementManager.h"
#include "../../../../analysis/FreqMNAM.h"
#include "NPort.h"

#include <cstdio>
#include <cstdlib>

extern "C"
{
  #include "../../../../compat/in_interp.h"
}
#include "../../../lex.NPort_.c"

// Static members
const unsigned NPort::n_par = 2;

// Element information
ItemInfo NPort::einfo =
{
  "nport",
  "Multi-port element defined by its port-based y parameters",
  "Carlos E. Christoffersen, Mark Summers",
  DEFAULT_ADDRESS"category:multiport",
  "2000_07_20"
};

// Parameter information
ParmInfo NPort::pinfo[] =
{
  {"filename",
	"File containing the port-based parameter matrix.", TR_STRING, true},
  {"max_freq",
	"Maximum number of frequency points in data file", TR_INT, false}
};


NPort::NPort(const string& iname) : Element(&einfo, pinfo, n_par, iname)
{
  // Standard construction code.
  paramvalue[0] = &filename;
  paramvalue[1] = &(max_freq = 200);

  setFlags(LINEAR | MULTI_REF | TR_FREQ_DOMAIN);

  // Set null pointers until we allocate memory in init()
  ymat_vec = NULL;
  freq_vec = NULL;

  // Set initial terminal and reference count from file.
  ports = 0;
  current_group = 0;
  // Set last_index to zero initially
  last_index = 0;
}


NPort::~NPort()
{
  if(ymat_vec)
	{
    // Free allocated matrices
    for (unsigned i=0; i < vec_size; i++)
      delete ymat_vec[i];

    delete [] ymat_vec;
    delete [] freq_vec;
  }
}

void NPort::check() throw(string&)
{
  // We need to check parameters first, before attepmting to open
  // the data file.
  try
	{
    checkParams();
  }
  catch(string& msg)
	{
    throw(getInstanceName() + ": " + msg);
  }
  if (max_freq < 1)
    throw(getInstanceName() + ": max_freq can not be less than 1");
}

void NPort::init() throw(string&)
{
  // Set initial dimension for internal frequency vector.
  vec_size = 0;
  // Now allocate memory for the vectors.
  freq_vec = new double[max_freq];
  ymat_vec = new ComplexDenseMatrix*[max_freq];

  FILE* fp;

  // Open input file
  if( (fp = fopen(filename.c_str(), "r")) == NULL)
    throw(getInstanceName() + ": Cannot open input file \"" +
	filename +"\".");

  // Point scanner pointer here.
  scan_NPort.thisNPort = this;
  NPort_in = fp;
  // Set the line counter to one.
  scan_NPort.linecount = 1;

  // Parse first part (ports) of input file
  if(NPort_lex())
    throw(getInstanceName() + ": in file \"" + filename + "\": " +
	+ scan_NPort.serror);

  // Check that the number of terminals in netlist is consistent
  // with number of terminals in file
  if (getTermCount() != ports + group_vec.size())
    throw(string(getInstanceName() + ": Incorrect number of terminals."));
  // If control reaches here, then the number of terminals in
  // netlist is correct.
  setNumTerms(getTermCount());

  // Resize the size of the workspace matrix
  ws_matrix.reshape(getNumTerms(), getNumTerms());

  // Continue parsing the rest of the file (y parameters).
  if (NPort_lex())
    throw(getInstanceName() + ": in file \"" + filename + "\": " +
	+ scan_NPort.serror);

  // Check that there is at least one set of parameters
  if (!vec_size)
    throw(getInstanceName() + ": in file \"" + filename +
	"\": At least one frequency data is required.");

  // Close input file
  fclose(fp);
}

void NPort::addPort(const unsigned& group_number)
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


ComplexDenseMatrix* NPort::newMatrix(const double& frequency)
{
  if (int(vec_size) == max_freq) {
    return NULL;
  }
  freq_vec[vec_size] = frequency;
  ymat_vec[vec_size] = new ComplexDenseMatrix(ports, ports);
  vec_size++;
  return ymat_vec[vec_size - 1];
}

void NPort::getLocalRefIdx(UnsignedVector& local_ref_vec,
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

void NPort::fillMNAM(FreqMNAM* mnam)
{
  const double& freq(mnam->getFreq());
  // Clean the workspace matrix.
  for (int i = 0; i < ws_matrix.numRows(); i++)
  {
    for (int j = 0; j < ws_matrix.numCols(); j++)
    {
      ws_matrix(i,j) = double_complex(zero,zero);
    }
  }
  // Find the frequency index.
  if (vec_size == 1)
    last_index = 0;
  else
    doubleHunt(freq_vec, int(vec_size), freq, &last_index);

  // First set the port submatrix into ws_matrix
  if (unsigned(last_index) >= vec_size)
	{
    // Desired frequency is too high
    // Use last information available.
    // We should also generate a warning message.
    last_index = vec_size - 1;

    for (unsigned i = 0; i < ports; i++)
    {
      for (unsigned j = 0; j < ports; j++)
      {
        ws_matrix(i,j) = (*ymat_vec[last_index])(i,j);
      }
    }
  }
  else if (last_index < 1)
	{
    // Use first matrix
    for (unsigned i = 0; i < ports; i++)
    {
      for (unsigned j = 0; j < ports; j++)
      {
        ws_matrix(i,j) = (*ymat_vec[0])(i,j);
      }
    }
  }
  else if (freq == freq_vec[last_index])
    // Desired frequency equals one of the data frequencies
	// Set element in MNAM
	for (unsigned i = 0; i < ports; i++)
  {
		for (unsigned j = 0; j < ports; j++)
    {
      ws_matrix(i,j) = (*ymat_vec[last_index])(i,j);
    }
  }
  else
    // Interpolate between last_index and last_index-1
	for (unsigned i = 0; i < ports; i++)
		for (unsigned j = 0; j < ports; j++)
		{
			// Calculate y parameter using interpolation
			double_complex yval((*ymat_vec[last_index])(i,j));
			yval -= (*ymat_vec[last_index-1])(i,j);
			yval *= ((freq - freq_vec[last_index-1]) /
			(freq_vec[last_index] - freq_vec[last_index-1]));
			yval += (*ymat_vec[last_index-1])(i,j);
			// Set element
      ws_matrix(i,j) = yval;
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

