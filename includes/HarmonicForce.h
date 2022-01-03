#pragma once



#include "output.h"
#include "PeriodicBoundary.h"
#include "MatMulVec.h"




class HarmonicForce{

	private:
		Real PotPar ;   // coupling constant
		std::vector<Real> EquiDist; // equilibrium distance

	public:
		HarmonicForce( Real x  , std::vector<Real> Delta ){ PotPar = x; EquiDist = Delta; };
		~HarmonicForce( void ){};
		std::vector<std::vector<Real>> ComputeForce( std::vector<std::vector<Real>> Pos , 
				   std::vector<std::vector<int>> LinkList );


		Real ComputePotential( std::vector<std::vector<Real>> Pos ,
                                 std::vector<std::vector<int>> LinkList  );






};

std::vector<std::vector<Real>> HarmonicForce::ComputeForce( std::vector<std::vector<Real>> Pos , std::vector<std::vector<int>> LinkList ){

	std::vector<std::vector<Real>> Force;
	Force.resize( Pos.size() );
	for ( auto i = 0 ; i < Pos.size() ; i++ ){
		Force[i].resize( Pos[i].size() );
		for ( auto j = 0 ; j < LinkList[i].size() ; j++ ){
			std::vector<Real> r = 
		     	        get_nearest_image( Pos[i] , Pos[LinkList[i][j]] ) - Pos[i];
			Real norm = Real( 1 ) / VecNorm( r );
			std::vector<Real> UnitVec = r * norm;
			std::vector<Real> delta = VecAbs( r ) - EquiDist;
			Force[i] = Force[i] - Real( 2 ) * PotPar * delta * UnitVec;
	        }
        }
	return Force;
 }

 Real HarmonicForce::ComputePotential( std::vector<std::vector<Real>> Pos ,
                                       std::vector<std::vector<int>> LinkList ){
	
	Real PotentialE = Real( 0.0 );
	for ( auto i = 0 ; i < Pos.size() ; i++ ){
		for ( auto j = 0 ; j < LinkList[i].size() ; j++ ){
			std::vector<Real> r = 
				get_nearest_image( Pos[i] , Pos[LinkList[i][j]] );
			std::vector<Real> delta = r - EquiDist;
			Real norm = VecNorm( delta );
			PotentialE = PotentialE  + PotPar * DotProduct( delta , delta );
	        }
        }
	return PotentialE;
 }

