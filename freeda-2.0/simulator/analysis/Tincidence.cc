// This function calculates T in compressed form and also returns the
// element vector, the number of state variables and the maximun
// number of state variables that an element can have.

#include "../network/Circuit.h"

void buildTIncidence(ElemFlag mask, Circuit*& my_circuit,
IntDenseMatrix& T, ElementVector& elem_vec,
int& n_states, int& max_n_states)
{
  n_states = 0;
  max_n_states = 0;
  // Clean vector
  elem_vec.clear();
  // Reserve a reasonable amount of memory.
  elem_vec.reserve(my_circuit->getNumberOfElements()/2);
  // Loop throw all selected elements in circuit.
  my_circuit->setFirstElement(mask);
  Element* elem = my_circuit->nextElement();
  while(elem)
  {
    int current_ns;
    // Ask the number of states aported by the element
    n_states += (current_ns = elem->getNumberOfStates());
    // Find the maximum number of states
    if (current_ns > max_n_states)
      max_n_states = current_ns;
    // Put the element in the vector
    elem_vec.push_back(elem);
    // get next element pointer
    elem = my_circuit->nextElement();
  }
  int n_elem = elem_vec.size();
  // Create T (incidence matrix).
  // T is frequency-independent, but the routine requires the row and
  // column numbers from the MNAM.
  // The format of T is a 2 x n_states matrix containing the row
  // indices for the 1 (row 0) and -1 elements (row 1).
  T.reshape(2, n_states);
  for (int m = 0; m < 2; m++)
    for (int n = 0; n < n_states; n++)
      T(m,n) = 0;
  int i = 0;
  for (int k = 0; k < n_elem; k++)
  {
    UnsignedVector local_ref_vec;
    TerminalVector term_list;
    // Get local reference node information
    elem_vec[k]->getLocalRefIdx(local_ref_vec, term_list);
    // This is just to make sure that things are consistent.
    assert(term_list.size() == elem_vec[k]->getNumTerms());
    // Number of terminal groups (local reference nodes).
    int ngroups = local_ref_vec.size();
    // jbase is the first terminal index in each group
    int jbase = 0;
    for (int l = 0; l < ngroups; l++)
		{
      int refcolumn = term_list[local_ref_vec[l]]->getRC();
      for (unsigned j = jbase; j < local_ref_vec[l]; j++)
			{
				// Add two matrix elements for each terminal: One element is a one
				// that goes in the terminal column, and the other is a minus one
				// in the local reference terminal column.
				int column = term_list[j]->getRC();
				// Remember that MNAM index begin at 1. Zero is ground.
				T(0,i) = column;
				T(1,i) = refcolumn;
				// change state, so we go to the next row.
				i++;
      }
      jbase = local_ref_vec[l] + 1;
    }
  }
}

