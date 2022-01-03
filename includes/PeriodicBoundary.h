#include <iostream>
#include <stdio.h>
#include <math.h>


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






#ifndef _PeriodicBound
#define _PeriodicBound

// returns sign of a val
template <typename T> int Sign(T val) {
    return ( T(0) < val ) - ( val < T(0) );
}




//
//returns nearest image of vector B with respect to A
std::vector<Real> get_nearest_image( std::vector<Real> A , std::vector<Real> B );
std::vector<Real> get_nearest_image( VecNSlice<Real> A , VecNSlice<Real> B );
std::vector<Real> get_nearest_image( std::vector<Real> A , VecNSlice<Real> B );
std::vector<Real> get_nearest_image( VecNSlice<Real> A , std::vector<Real> B );

//
// return nearest image of vector B with respect to vector A
// in cartesian coordinates
//
std::vector<Real> get_nearest_image_cart( std::vector<Real> A , std::vector<Real> B ,
		           std::vector<std::vector<Real>> lattice );
void ApplyPeriodicBound( std::vector<Real>& Pos );
void ApplyPeriodicBound( VecNSlice<Real>& Pos );

#endif
