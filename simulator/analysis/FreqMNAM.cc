#include "FreqMNAM.h"
#include <iostream>
using std::cout;
using std::endl;
using std::cerr;


// ------------------------------------------------------------
// Constructor
// ------------------------------------------------------------
FreqMNAM::FreqMNAM(const DenseDoubleVector& fvec, Circuit* cir,
const ElemFlag& mask)
{
  // resize freq_vec to the size of the fvec so
  // as to prevent segmentation faults
  freq_vec.resize(fvec.length());
  for (unsigned int i = 0; i < fvec.length(); i++)
    freq_vec[i] = fvec[i];
  current_findex = 0;

  // Set internal variables first.
  assert(cir);
  my_circuit = cir;
  mymask = mask;

  // Fill the element vector.
  // Reserve a reasonable amount of memory.
  // This routine allocates some default memory to hold the elements
  // in the Element vector
  elem_vec.reserve(my_circuit->getNumberOfElements());

  // Loop throw all selected elements in circuit.
  // the circuit pointer is pointing to the first element
  my_circuit->setFirstElement(mask);

  // ..and is now looping through all the elements in the circuit
  // if it finds new elements during looping, it adds it to the element vector
  Element* elem = my_circuit->nextElement();
  while(elem)
  {
    // Put the element in the vector
    elem_vec.push_back(elem);
    // get next element pointer
    elem = my_circuit->nextElement();
  }

  // Local variables
  // Number of elements, the size of the element vector just created
  unsigned n_elem = elem_vec.size();

  // The number of frequencies, the size of the frequency vector
	// which is a constructor argument
  unsigned n_freqs = freq_vec.length();

  // Scratch variables
  unsigned findex;
  int error;

  // Set the initial value for the MNAM dimension.
  dimension = 0;

  // Set row/column information on terminals
  // Also find the reference node(s)
  // Now it is time to expand the MNAM depending on the no of terminals
	// present. The circuit pointer goes through all the terminals and
	// simultaneously searches for the reference terminal. If a terminal
	// found is not a reference temrinal then the dimension of the matrix is
	// increased. If the terminal is a reference terminal, then the boolean
	// variable "found" is set to true.
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
  // the getExtraRC function is defined in Element.cc and it's function is to
  // accomodate any extra equations that the element may require. It finds out
  // this extra number and returns this number which is then added to our
  // dimension varable, so that the correct matrix size is created
  for (unsigned i=0; i < n_elem; i++)
    dimension += elem_vec[i]->getExtraRC(dimension + 1, FREQ_DOMAIN);

  // Allocate space for the MNAM matrices, one for each frequency
  SerialCommunicator FreqMNAMComm;
  int nnz = (int)dimension/4.0;
  ProcessorMap Map(dimension, 0, FreqMNAMComm);

  for (int i = 0; i < n_freqs; i++)
  {
    // dynamically create the sparse A[i]
    Teuchos::RCP< DoubleSparseColMatrix > matReal = 
        Teuchos::rcp(new DoubleSparseColMatrix(Copy, Map, nnz));
    Teuchos::RCP< DoubleSparseColMatrix > matImag = 
        Teuchos::rcp(new DoubleSparseColMatrix(Copy, Map, nnz));
    // add it to the associated vectors
    Areal.push_back(matReal);
    Aimag.push_back(matImag);
  }

  // Allocate memory for the vector of source vectors.
  // For this we need to know first the real dimension of the MNAM
  // Initialize the vectors to zero.
  // The source vector size depends on the number of frequencies specified
  // and is created accordingly. For each frequency there must be a 'new'
  // source vector and whose dimension will be equal to the dimension variable
  // calculated above.
  sourceV_vec = new DenseComplexVector*[n_freqs];
  for (findex = 0; findex < n_freqs; findex++)
    sourceV_vec[findex] = new DenseComplexVector(dimension);
  // The source vector is initialized automatically
  // to zero at each frequency

  // Now fill the matrices, tell each element to fill the MNAM
  // At each frequency, the MNAM is filled with the value in the Element vector
  for (current_findex = 0; current_findex < n_freqs; current_findex++)
    for (unsigned i=0; i < n_elem; i++)
      elem_vec[i]->fillMNAM(this);
}


// ------------------------------------------------------------
// Destructor
// ------------------------------------------------------------
FreqMNAM::~FreqMNAM()
{
  for (unsigned findex = 0; findex < freq_vec.length(); findex++)
    delete sourceV_vec[findex];
  delete [] sourceV_vec;
}


// ------------------------------------------------------------
// explicitly declare the MNAMs complete, and build the augmented matrix
// ------------------------------------------------------------
void FreqMNAM::FillComplete()
{
  // now that all the Areal's and Aimag's are filled, call their respective
  // FillComplete() functions.
  unsigned n_freqs = freq_vec.length();
  for (int i = 0; i < n_freqs; i++)
  {
    (*Areal[i]).FillComplete();
    (*Aimag[i]).FillComplete();
  }
  FreqMNAM::buildA();
}


// ------------------------------------------------------------
// build the augmented matrix A
// ------------------------------------------------------------
void FreqMNAM::buildA()
{
  SerialCommunicator FreqMNAMComm2;
  unsigned n_freqs = freq_vec.length();
  int nnz = (int)dimension/4.0;
  ProcessorMap Map2(dimension*2, 0, FreqMNAMComm2);
  
  for (int i = 0; i < n_freqs; i++)
  {
    // dynamically create the sparse A[i]
    Teuchos::RCP< DoubleSparseColMatrix > matA = 
        Teuchos::rcp(new DoubleSparseColMatrix(Copy, Map2, nnz*4));
    // add it to the associated vectors
    A.push_back(matA);
  }
  
  // now build augmented matrix A
  double * rowValExtract = new double[dimension];
  int * colIndExtract = new int[dimension];
  int num_non_zeros;

  for (int i = 0; i < n_freqs; i++)
  {
    for (int k = 0; k < dimension; k++)
    {
      (*Areal[i]).ExtractGlobalRowCopy(k, dimension, num_non_zeros, 
        rowValExtract, colIndExtract);
      (*A[i]).InsertGlobalValues(k, num_non_zeros, rowValExtract, colIndExtract);
      for (int j = 0; j < num_non_zeros; j++)
        colIndExtract[j] += dimension;
      (*A[i]).InsertGlobalValues(k + dimension, num_non_zeros, 
        rowValExtract, colIndExtract);
    }

    for (int k = 0; k < dimension; k++)
    {
      (*Aimag[i]).ExtractGlobalRowCopy(k, dimension, num_non_zeros, 
        rowValExtract, colIndExtract);
      (*A[i]).InsertGlobalValues(k + dimension, num_non_zeros, 
        rowValExtract, colIndExtract);
      for (int j = 0; j < num_non_zeros; j++)
      {
        colIndExtract[j] += dimension;
        rowValExtract[j] = -rowValExtract[j];
      }
      (*A[i]).InsertGlobalValues(k, num_non_zeros, rowValExtract, colIndExtract);
    }
    (*A[i]).FillComplete();
  }
  delete [] rowValExtract;
  delete [] colIndExtract;
}


// ------------------------------------------------------------
// Methods related to the frequency vector
// ------------------------------------------------------------
void FreqMNAM::setFreqIndex(const unsigned& findex)
{
  assert(findex < freq_vec.length());
  current_findex = findex;
}


// ------------------------------------------------------------
// Methods to fill the matrices used with the elements
// ------------------------------------------------------------
// add a complex element to the MNAM at position (row,col)
void FreqMNAM::setElement(const unsigned& row, const unsigned& col, 
  const double_complex& val)
{
  double * rowVal = new double[1];
  int * colInd = new int[1];
  // check that row and col are non-negative and not out of bounds
  if ((row < 0) || (col < 0) || (row > dimension) || (col > dimension))
  {
    cerr << "Incorrect row or column position in FreqMNAM requested." << endl;
    exit(1);
  }
  
  // if row or col number is zero, then add nothing to the MNAM
  if ((row == 0) || (col == 0))
    return;
  
  // Now that all checks are complete, add to the MNAM.
  // The complex number is split into the real and imaginary part.
  // The real part is added to Areal, the imaginary part to Aimag
  // This is so that the augmented matrix A, created in the constructor above,
  // consists of real numbers only although it is of twice the dimension.
  // Currently, there is no generic support for complex numbers in Trilinos
  // and so this augmented real-number form must be created.
  rowVal[0] = val.real();
  colInd[0] = col-1;
  (*Areal[current_findex]).InsertGlobalValues(row-1, 1, rowVal, colInd);
  rowVal[0] = val.imag();
  (*Aimag[current_findex]).InsertGlobalValues(row-1, 1, rowVal, colInd);
  
  delete [] rowVal;
  delete [] colInd;
}


// ------------------------------------------------------------
// this function adds the admittance value to four locations in the MNAM
// and for the addition at each location, it calls the setElement() above.
// At location (term1, term1) and (term2, term2) it adds +val
// At location (term1, term2) and (term2, term1) it adds -val
// ------------------------------------------------------------
void FreqMNAM::setAdmittance(const unsigned& term1, const unsigned& term2,
  const double_complex& val)
{
  // no checks need be performed here, since setElement() above will do it
  FreqMNAM::setElement(term1, term1, val);
  FreqMNAM::setElement(term2, term2, val);
  FreqMNAM::setElement(term1, term2, -val);
  FreqMNAM::setElement(term2, term1, -val);
}


// ------------------------------------------------------------
// this function is a more general form of setAdmittance() where the row
// and col can be any numbers within the range of the MNAM dimensions, 
// and do not have to be a pair. For each location, setElement() is called
// At location (row1, col1) it adds +val
// At location (row2, col2) it adds +val
// At location (row1, col2) it adds -val
// At location (row2, col1) it adds -val
// ------------------------------------------------------------
void FreqMNAM::setQuad(const unsigned& row1, const unsigned& row2,
  const unsigned& col1, const unsigned& col2, const double_complex& val)
{
  // no checks need be performed here, since setElement() above will do it
  FreqMNAM::setElement(row1, col1, val);
  FreqMNAM::setElement(row2, col2, val);
  FreqMNAM::setElement(row1, col2, -val);
  FreqMNAM::setElement(row2, col1, -val);
}


// ------------------------------------------------------------
// this function adds a positive and negative one to the MNAM at locations
// indicated by eqn. It makes calls to setElement() to accomplish this.
// At location (neg, eqn) it adds a -1
// At location (eqn, neg) it adds a -1
// At location (pos, eqn) it adds a +1
// At location (eqn, pos) it adds a +1
// ------------------------------------------------------------
void FreqMNAM::setOnes(const unsigned& pos, const unsigned& neg,
  const unsigned& eqn)
{
  // some additional checks need to be performed here before calling
  // setElement()
  // eqn cannot be 0. If it is, then add nothing to the MNAM
  if (eqn == 0)
    return;
  // pos and neg cannot be of the same value
  if (pos == neg)
    return;
  // neither pos or neg should be the same as eqn
  if ((pos == eqn) || (neg == eqn))
    return;
  
  FreqMNAM::setElement(neg, eqn, double_complex(-1.0));
  FreqMNAM::setElement(eqn, neg, double_complex(-1.0));
  FreqMNAM::setElement(pos, eqn, double_complex(1.0));
  FreqMNAM::setElement(eqn, pos, double_complex(1.0));
}


// ------------------------------------------------------------
// It gets the source vector at a particular frequency
// ------------------------------------------------------------
void FreqMNAM::getSource(const unsigned& findex, DenseComplexVector& sourceV)
{
  assert(findex < freq_vec.length());
	if (sourceV.length() != dimension)
		sourceV.resize(dimension);
  // copy sourceV_vec[findex] to source
  for (int i = 0; i < dimension; i++)
    sourceV[i] = (*sourceV_vec[findex])[i];
}


// ------------------------------------------------------------
// Factors and Solve a matrix system for a fixed source vector, i.e. for a
// particular frequency. The solution is stored in the result vector.
// Here the size of the result vector is set to the value 'dimension'
// ------------------------------------------------------------
void FreqMNAM::solve(const unsigned& findex, DenseComplexVector& result)
{
  assert(findex < freq_vec.length());
  // Check dimension of result vector
  if (result.length() != dimension)
    result.resize(dimension);
  
  SerialCommunicator Comm2;
  ProcessorMap Map2(dimension*2, 0, Comm2);
  
  // X and B must also be twice the size, because they contain 
  // real components only
  DistributedDoubleVector X(Map2);
  DistributedDoubleVector B(Map2);
  
  for (int i = 0; i < dimension; i++)
  {
    B[i] = (*sourceV_vec[findex])[i].real();
    B[i+dimension] = (*sourceV_vec[findex])[i].imag();
  }

  // have to factorize first before solving
  Problem.SetOperator(&*A[findex]);
  Problem.SetLHS(&X);
  Problem.SetRHS(&B);
  const char * SolverType = "Klu";
  Solver = Factory.Create(SolverType, Problem);
  // perform factorization ...
  int ierr = Solver->NumericFactorization();
  if (ierr > 0)
  {
    cerr << "**Problem factoring the MNAM!**" << endl;
    exit(1);
  }
    
  // ... and back substitution to put solution into X
  Solver->Solve();
  
  // Finally, copy the solution from X into result
  for (int i = 0; i < dimension; i++)
  {
    result[i].real() = X[i];
    result[i].imag() = X[i+dimension];
  }
}


// ------------------------------------------------------------
// It solves a matrix system for a given source vector
// The solution is stored in the result vector.
// Here the size of the result vector is set to the value 'dimension'
// ------------------------------------------------------------
void FreqMNAM::solve(const unsigned& findex, const DenseComplexVector& sourceV,
DenseComplexVector& result)
{
	assert(findex < freq_vec.length());
  assert(sourceV.length() == dimension);
  // Check dimension of result vector
  if (result.length() != dimension)
    result.resize(dimension);
  
  SerialCommunicator Comm2;
  ProcessorMap Map2(dimension*2, 0, Comm2);
  
  // X and B must also be twice the size, because they contain 
  // real components only
  DistributedDoubleVector X(Map2);
  DistributedDoubleVector B(Map2);
  
  for (int i = 0; i < dimension; i++)
  {
    B[i] = sourceV[i].real();
    B[i+dimension] = sourceV[i].imag();
  }

  // have to factorize first before solving
  Problem.SetOperator(&*A[findex]);
  Problem.SetLHS(&X);
  Problem.SetRHS(&B);
  const char * SolverType = "Klu";
  Solver = Factory.Create(SolverType, Problem);
  // perform factorization ...
  int ierr = Solver->NumericFactorization();
  if (ierr > 0)
  {
    cerr << "**Problem factoring the MNAM!**" << endl;
    exit(1);
  }

  // ... and back substitution to put solution into X
  Solver->Solve();
  
  // Finally, copy the solution from X into result
  for (int i = 0; i < dimension; i++)
  {
    result[i].real() = X[i];
    result[i].imag() = X[i+dimension];
  }
}
  

// ------------------------------------------------------------
// Method to fill the source vectors
// ------------------------------------------------------------
void FreqMNAM::addToSource(const unsigned& pos, const unsigned& neg,
const double_complex& val)
{
  assert(pos <= dimension);
  assert(neg <= dimension);
  assert(current_findex < freq_vec.length());
  if (pos)
    (*sourceV_vec[current_findex])[pos-1] += val;
  if (neg)
    (*sourceV_vec[current_findex])[neg-1] -= val;
}

