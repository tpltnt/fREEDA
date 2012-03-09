#include <cassert>
#include "TimeMNAM.h"

TimeMNAM::TimeMNAM(Circuit* cir, const ElemFlag& mask) :
curr_time(0), mymask(mask), my_circuit(cir)
{
  // Set internal variables first.
  assert(cir);
	
  // Fill the element vector.
  // Reserve a reasonable amount of memory.
  elem_vec.reserve(my_circuit->getNumberOfElements());
  // Loop throw all selected elements in circuit.
  my_circuit->setFirstElement(mask);
  Element* elem = my_circuit->nextElement();
  while(elem) 
	{
		// Put the element in the vector
		elem_vec.push_back(elem);
		// get next element pointer
		elem = my_circuit->nextElement();
	}
	
  // Now fill the source element vector
  // Reserve a reasonable amount of memory.
  source_elem_vec.reserve(elem_vec.size());
  // Loop throw all selected source elements in circuit.
  my_circuit->setFirstElement(mask | SOURCE);
  elem = my_circuit->nextElement();
  while(elem) 
	{
		// Put the element in the vector
		source_elem_vec.push_back(elem);
		// get next element pointer
		elem = my_circuit->nextElement();
	}
	
  // Local variables
  // Number of elements
  unsigned n_elem = elem_vec.size();
	
  // Set the initial value for the MNAM dimension.
  dimension = 0;
	
  // Set row/column information on terminals
  // Also find the reference node(s)
  my_circuit->setFirstTerminal();
  Terminal* term = my_circuit->nextTerminal();
  bool found = false;
  while (term) 
	{
		// check if the terminal is a reference
		if (term->isRef()) 
		{
			term->setRC(0);
			found = true;
		}
		else 
		{
			// The dimension is the non-reference terminal count.
			dimension++;
			term->setRC(dimension);
		}
		// get next terminal
		term = my_circuit->nextTerminal();
	}
  assert(found);
	
  // Get final dimension asking for extra rows and columns
  for (unsigned i=0; i < n_elem; i++) 
    dimension += elem_vec[i]->getExtraRC(dimension + 1, TIME_DOMAIN);
	
  // create matrices and source vector
  int dim(dimension);
  SerialCommunicator TimeMNAMComm;
  ProcessorMap Map(dim, 0, TimeMNAMComm);
  M = new DoubleSparseColMatrix(Copy, Map, int(dim/2));
  Mp = new DoubleSparseColMatrix(Copy, Map, int(dim/2));
  // the sizes of rValues and colIndices are set to one because we want to
  // set individual values in the sparse matrices
  rowVal = new double[1];
  colInd = new int[1];
  Lhead = new_sllist();
  clist = Lhead;
  nlist = Lhead->next;
  // Now fill the matrices
  for (unsigned i=0; i < n_elem; i++)
  {
    // Tell each element to fill the MNAM
	  elem_vec[i]->fillMNAM(this);
  }
}

TimeMNAM::~TimeMNAM()
{
  delete[] rowVal;
  delete[] colInd;
  delete Mp;
  delete M;
}

void TimeMNAM::getSource(const double& time, DenseDoubleVector& ansV)
{
  // ask the sources to fill the vector
  curr_time = time;
  // Point our source vector to the
  sourceV = &ansV;
  // Clear vector
  sourceV->putScalar(zero);
  unsigned n_elem = source_elem_vec.size();
  //assert(ansV.length() >= n_elem);
	// Tell each element to fill the MNAM
  for (unsigned i=0; i < n_elem; i++)
    source_elem_vec[i]->fillSourceV(this);
}

void TimeMNAM::print() const
{
  M->Print(std::cout);
  cout << "\n\n\n";
  Mp->Print(std::cout);
  //M1->Print(std::cout);
}
