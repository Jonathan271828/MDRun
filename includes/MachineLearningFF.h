#pragma once
#include <vector>
#include <fstream>
#include <string>
#include <iostream>
#include <string>
#include <sstream>
#include <numeric>
#include <vector>
#include <ios>
#include <functional>
#include <math.h>



#ifndef _VectorTools
#include "MatMulVec.h"
#endif


#ifndef _FLAGFinder
#include "ReadFlags.h"
#endif


#ifndef _MDSPEC
#include "atoms.h"
#endif



#ifndef _PeriodicBound
#include "PeriodicBoundary.h"
#endif

#ifndef _LinkList
#include "LinkList.h"
#endif





#ifdef USE_DOUBLES
typedef double Real;
#else
typedef float Real;
#endif





/*
 *
 *
 * object to read tables from file
 * first call constructor with 
 * filename you want to read from
 * call ReadData
 * and last call ReturnData
 * to get read data as 2d array
 *
 *
 *
 */


template <class T>
class ReadTable{

	private:
		std::ifstream infile;
		std::vector<std::vector<T> > data;

		std::vector<T> ExtractValues( std::string input ){
        
        	std::string Temp;
        	std::vector<T> Data ;
        
        	std::stringstream StrStream( trim( input ).c_str() );
        
        
        	while ( StrStream >> Temp ){
        		Data.push_back( ::atof( Temp.c_str() ) );
        	}
			return Data;
		}


	    

	public:
        ReadTable( std::string fname ){
			infile.open ( fname.c_str(), std::ifstream::in );
        };

        void ReadData( void ){
        
        	std::string Filedata;
        	std::string x;
        
        	infile.clear();
        	infile.seekg( 0 , std::ios::beg );
        
        	while ( std::getline( infile , Filedata ) ){
				data.push_back( ExtractValues( Filedata ) );
        	}
        
        };

		std::vector<std::vector<T> > ReturnData( void ){
			return data;
		}


};







class MachineLearningForces{

	private:
		std::string WeightFile;
		std::string DescriptorFile;
		std::string DescriptorDerivativeFile;

		std::vector<Real> Weights;
	    std::vector<VecNSlice<Real> > Descriptors;
		std::vector<VecNSlice<Real> > DescriptorsDerivative;
		std::vector<VecNSlice<Real> > DescriptorActStruc;
		std::vector<VecNSlice<Real> > DescriptorDerivativeActStruc;
	    std::vector<Real> DescriptorsData;
		std::vector<Real> DescriptorsDerivativeData;
		std::vector<Real> DescriptorActStrucData;
		std::vector<Real> DescriptorDerivativeActStrucData;

		std::vector<Real> Frequencies;
		std::vector<Real> FilterFunction;



		unsigned int NBasis ;


	public:
		ComputeLinkLists LinkList;

	    MachineLearningForces( std::string file1 , std::string file2 , std::string file3 , 
				               std::function <void ( std::vector<Real>& data , std::vector<VecNSlice<Real> >& Pointers ,
				   	           unsigned int dim )> y , unsigned int Natoms , Real Cutoff ,  
							   std::vector<VecNSlice<Real> > PosDir , std::vector<VecNSlice<Real> > PosCar ,
							   std::vector<std::vector<Real> > Box );
       void ComputeDescriptors( std::vector<VecNSlice<Real> >& PosDir , const std::vector<std::vector<Real> >& Box );
       
	   void ComputeForces( std::vector<VecNSlice<Real> >& ForceDir , std::vector<VecNSlice<Real> >& ForceCar ,
			               const std::vector<std::vector<Real> > Box );
	   Real ComputeEnergy( void );

};
