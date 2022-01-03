#pragma once


#include "PeriodicBoundary.h"
#include "output.h"



#ifdef USE_DOUBLES
typedef double Real;
#else
typedef float Real;
#endif





 class MexicanHat
 {

	public:
		Real PotParA; //parameters for potential
		Real PotParB; //parameters for potential
		std::vector<std::vector<Real>> CoP;  // center of potential


		MexicanHat( Real A = Real( 1 ) , Real B = Real( 1 ) , 
	         std::vector<std::vector<Real>> Centers = 
		 std::vector<std::vector<Real>> ( 3 , std::vector<Real>( 3 , Real( 0 ) )  ) );
		~MexicanHat( void ){};
		Real ComputePotentialEnergy( std::vector<std::vector<Real>> Pos , std::vector<Real> Epot );
		void ComputeForces( std::vector<std::vector<Real>> Pos , std::vector<std::vector<Real>>& Force );

 };


 MexicanHat::MexicanHat( Real A , Real B , std::vector<std::vector<Real>> Centers ){
	 PotParA = A ;
	 PotParB = B ; 
 }


 Real MexicanHat::ComputePotentialEnergy( std::vector<std::vector<Real>> Pos , std::vector<Real> Epot ){


	 Real Etot = Real( 0 );


	 for ( auto i = 0 ; i < Pos.size() ; i++ ){
		 std::vector<Real> dr = SubVecs( CoP[i] , 
				 get_nearest_image( CoP[i] , Pos[i] ) );
		 Real norm = VecNorm( dr );
		 Epot[ i ]  =  PotParA*norm*norm  +  PotParB*norm*norm*norm*norm;
		 Etot = Etot + Epot[i];
	 }

	 return Etot;
 }

 
 void MexicanHat::ComputeForces( std::vector<std::vector<Real>> Pos , std::vector<std::vector<Real>>& Force ){


	 /// minus for force is taken into account during SubVecs function call
	 for ( auto i = 0 ; i < Pos.size() ; i++ ){
		 std::vector<Real> dr = SubVecs( CoP[i] , get_nearest_image( CoP[i] , Pos[i] ) );
		 Real norm  =  VecNorm( dr );
		 for ( auto j = 0 ; j < Pos[i].size() ; j++ ){
			 Force[i][j] =  Real( 2 )*PotParA*dr[j] + 
				        Real( 4 )*PotParB*norm*norm*dr[j];
		 }
	 }
 }



