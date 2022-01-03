#include <iostream>
#include <math.h>
#include <vector>
#include <fstream>
#include <string>
#include <sstream>
#include <numeric>
#include <ctime>
#include <complex>
#include <Eigen/Eigenvalues>
#include <iomanip>
#include <ios>

#ifndef _VectorTools
#include "MatMulVec.h"
#endif
#ifndef _Forces
#include "force.h"
#endif
#ifndef _FLAGFinder
#include "ReadFlags.h"
#endif
#ifndef _PeriodicBound
#include "PeriodicBoundary.h"
#endif

#ifndef _MDSPEC
#include "atoms.h"
#endif



#ifdef USE_DOUBLES
typedef double Real;
#else
typedef float Real;
#endif


typedef std::complex<Real> CReal;


std::vector<Real> extractVector( std::string data );

#ifndef _Phononsfidi
#define _Phononsfidi
  class FDPhonons{

	  public:
		  // initializer reads qpoints file and computes the unit cells; 
		  // expresses Qpoints supplied in cartesian coordinates
		  // in cartesian coordinates of the Reciprocal unit cell times a factor two pi
		  FDPhonons( std::vector<std::vector<Real>> Box , 
				  std::vector<std::vector<Real>> RecBox , std::vector<int> NN ,
				  std::vector<Real> Dx );
		  // main calls all the necessary routines ; it needs positions in 
		  // direct and cartesian coordinates the class the computes the force, 
		  // the system Box, and the atomic masses
		  void main( std::vector<VecNSlice<Real> >& Pos , 
			     ComputeForces & GetForce ,
			     std::vector<std::vector<Real> >& Box ,
			     std::vector<Real>& Mass );


	  private:
		  // units needed for quantum espresso output
		  const Real eVtTHz = 241.79893;  // conversion from electron volts to terahertz
		  const Real eVtIcm = 8065.54429; // conversion from electron volts to cm^-1

		  std::vector<std::vector<Real> >  Qpoints;  // array containing Qpoints
		  std::vector<std::vector<Real> >  UnitCell;
		  std::vector<std::vector<Real> >  RecUnitCell;  // reciprocal unit cell (no factor 2pi included)
		  std::vector<Real> delta;   // defines the magnitude of the atomic displacements
		  std::vector<std::vector<std::vector<std::vector<std::vector<std::vector<Real> > > > > >
			                 ForceConstMatrix;
		  std::vector<Eigen::MatrixXcd> DynamicMatrix;
//		  Eigen::ComplexEigenSolver<Eigen::MatrixXcd> ComplexSolver;  // solver for eigenvalues and eigenvectors of dynamic matrix
                  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> ComplexSolver;
		  std::vector<std::vector<CReal> >EigenValues;
		  std::vector<std::vector<std::vector<CReal> > > EigenVectors;
		  std::ofstream DebugFile;
		  // compute the force constant matrix in real space by central differences
                  void ComputeForceConstantMatrix( std::vector<VecNSlice<Real> >& Pos ,
				  ComputeForces& GetForce , 
				  std::vector<std::vector<Real>>& Box , std::vector<Real>& Mass );
		  // computes dynamical matrix via fourier transform
		  void ComputeDynamicMatrix( std::vector<VecNSlice<Real> > Pos , std::vector<std::vector<Real> > Box );
		  // compute eigenvalues and eigenvectors with eigensolver and strore them to the arrays EigenValues and EigenVectors
		  void ComputeEigenVectorsValues( void );
		  // writing the computed data to a files
		  void WriteOutput( void );
		  // apply acoustic sum rule
		  void AcousticSumRule( unsigned int dim );
		  // symmetrize dynamic matrix after fourier transforming
		  void Symmetrize( unsigned int dim );
		  // compute second derivative of potential
		  Real ComputeSecondDerivativeOD( VecNSlice<Real> PosA , VecNSlice<Real> PosB , 
				                std::vector<std::vector<Real>> Box, 
						unsigned int alpha , unsigned int beta , ComputeForces GetForce );
		  // computes the second derivative for equal atoms
                  Real ComputeSecondDerivativeDE( std::vector<VecNSlice<Real> >& PosA , 
				                  std::vector<std::vector<Real> >& Box,
		                                  unsigned int alpha, unsigned int beta , 
						  ComputeForces & GetForce ,
					          unsigned int AtNr );
		  // write dyn.out file in quantum espresso format
		  // as input for the phq code
		  void WriteQuantumEspressoOutput( void );
		  // write dyn.out file in a very simplified way
		  void WriteSimpleOutput( void );
		  // write special entry of dynamic matrix
		  void PrintDynMatrix( unsigned int QN );
		  void WriteForceMatrix( void );
  };

#endif
