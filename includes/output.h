
#include <iostream>
#include <stdio.h>
#include <string>
#include <vector>



#ifndef _VectorTools
#include "MatMulVec.h"
#endif

#ifndef _MDSPEC
#include "atoms.h"
#endif



#ifdef USE_DOUBLES
typedef double Real;
#else
typedef float Real;
#endif



#ifndef _OUTPUT
#define _OUTPUT


//
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// ~~~~~~~~~~~~~~~~~~~   PRINT ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

 template <typename T>
 void Print( std::vector<T> out ){

	for ( auto i = 0 ; i < out.size() ; i++ ){
		std::cout << "   " << out[i];
	}
	std::cout << std::endl;
 }

 template <typename T>
 void Print( T x ){
 	std::cout << x << std::endl;
 }

 template <typename T>
 void Print( std::vector<std::vector<Real>> A ){

	 for ( auto i = 0 ; i < A.size() ; i++ ){
		 for ( auto j = 0 ; j< A[i].size(); j++ ){
		    std::cout << "   " << A[i][j];
		 }
	     std::cout << std::endl;
	 }
 }



 void ErrorMessage( std::string text , std::string routine );

// write print function here that is calling the others with variable
// argument number


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//







class WriteOutput{


	private:
	FILE *file;
	unsigned int Counter = 0;

	public:
	WriteOutput( FILE * fp ){ file = fp; };
	void XdatInit( std::vector<std::vector<Real>> lattice , std::vector<std::string> Spec ,
			 std::vector <int> Natoms );
	void AddStruc( std::vector<VecNSlice<Real>> Pos );
	void WriteVector( std::vector<Real> data );
	void LammpsWrite( std::vector<VecNSlice<Real> > Pos ,
			  std::vector<std::vector<Real>> lattice , std::vector<int> Spec ,
                               unsigned int Natoms );
	void PhQWriteStruc( std::vector<VecNSlice<Real> > Pos , bool first , unsigned int Number );
};

#endif
