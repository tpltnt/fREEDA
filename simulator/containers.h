// The interface file for containers used by fREEDA

#ifndef containers_h
#define containers_h 1

#include <cstring>
#include <complex>
#include <vector>
#include <cmath>

#ifdef HAVE_MPI
#include "mpi.h"
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif

#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_CrsMatrix.h"
#include "Teuchos_SerialDenseMatrix.hpp"
#include "Teuchos_SerialDenseVector.hpp"

typedef std::complex<double> double_complex;

// ----------------------------------------
// typedefs for all vector containers
// ----------------------------------------
// Dense vectors of type double, int and complex
typedef Teuchos::SerialDenseVector<int,double> DenseDoubleVector;
typedef Teuchos::SerialDenseVector<int,int> DenseIntVector;
typedef Teuchos::SerialDenseVector<int,double_complex> DenseComplexVector;

// Distributed vectors of type double
typedef Epetra_Vector DistributedDoubleVector;
typedef Epetra_MultiVector DistributedDoubleMultiVector;

// STL vector
typedef std::vector<unsigned> UnsignedVector;

// ----------------------------------------
// typedefs for all matrix containers
// ----------------------------------------
// Sparse column oriented matrix
typedef Epetra_CrsMatrix DoubleSparseColMatrix;
typedef Epetra_SerialComm SerialCommunicator;
typedef Epetra_Map ProcessorMap;

// Dense matrices of type double, int and complex 
typedef Teuchos::SerialDenseMatrix<int,double> DoubleDenseMatrix;
typedef Teuchos::SerialDenseMatrix<int,int> IntDenseMatrix;
typedef Teuchos::SerialDenseMatrix<int,double_complex> ComplexDenseMatrix;

// ----------------------------------------
// declare some handy global constants
// ----------------------------------------
extern const double zero;
extern const double one;
extern const double_complex Czero;
extern const double_complex Cone;
extern const double pi;
extern const double twopi;
// Constant to convert degrees to radians (divide for inverse conversion)
extern const double deg2rad;
// Boltzman constant
extern const double kBoltzman;
// Electron charge
extern const double eCharge;
// Permittivity of free space 
extern const double epsilon0;
// Permeability of free space
extern const double mu0;
// Speed of light in free space
extern const double c0;

#endif

