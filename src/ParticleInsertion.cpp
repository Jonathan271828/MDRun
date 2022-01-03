
#include "ParticleInsertion.h"






void InsertParticle::DefineLinkList( std::vector<Real> PosNew , std::vector<std::vector<Real> > Pos , Real Cutoff , 
		                             std::vector<std::vector<Real> > Box ){
	
	LinkList.resize( 0 );
	for ( auto i = 0 ; i < Pos.size() ; i++ ){
		Real norm = VecNorm( MatVec( Box , get_nearest_image( Pos[i] , PosNew ) - PosNew ) );
		if ( norm < Cutoff ){
			LinkList.push_back( i );
		}
	}
}


std::vector<Real> InsertParticle::MonteCarloInsertion( std::vector<VecNSlice<Real> > & Pos , 
		                                       std::vector<VecNSlice<Real> > & Velos , const Real Cutoff , 
		                                       const std::vector<std::vector<Real> > Box , \
						       ComputeForces force ,
						       Real Sigma ){
	

//	std::vector<Real> PosOldArray = NewPosition( Pos[0].size() );
////	VecNSlice<Real> PosOld =  PosOldArray;
//
//
//
//	DefineLinkList( PosOld , Pos , Cutoff , Box );
//	Real energyOld = force.SupplyMorseSingleAtom( PosOld , Pos , LinkList , Box ) ;
//
//	std::vector<Real> forceNew = force.SupplyMorseForceSingleAtom( PosOld , Pos , LinkList );
//	for ( auto i = 0 ; i < Nsteps ; i++ ){
//		std::vector<Real> PosNew  =  NewPosition( Pos[0].size() );
//		DefineLinkList( PosNew , Pos , Cutoff , Box );
//		Real energyNew = force.SupplyMorseSingleAtom( PosNew , Pos , LinkList , Box );
//		if ( energyNew < energyOld ){
//			PosOld = PosNew;
//			energyOld = energyNew;
//			forceNew = force.SupplyMorseForceSingleAtom( PosNew , Pos , LinkList );
//
//		}
//		else {
//			Real rand   =  Rand.doub();
//			Real deltaE =  energyNew - energyOld;
//			Real Boltz =  exp( -deltaE / temperature );
//			if ( rand < Boltz ){
//				energyOld = energyNew;
//				PosOld = PosNew;
//			    forceNew = force.SupplyMorseForceSingleAtom( PosNew , Pos , LinkList );
//			}
//		}
//	    std::cout << "Step " << i << "  " << energyOld << std::endl;
//	}
//	// velocity mak e
//	std::vector<Real> VeloNew;
//	for ( auto i = 0 ;  i < Pos[0].size() ; i++ ){
//		VeloNew.push_back( Rand.generate_gauss_dist( Real( 0 ) , Sigma ) );
//	}
//	Velos.push_back( VeloNew );
//	Pos.push_back( PosOld );
//	return forceNew;
}



std::vector<Real> InsertParticle::NewPosition( unsigned int N ){
	std::vector<Real> Pos;
	for ( auto i = 0 ; i < N ; i++ ){
		Pos.push_back( Rand.doub() );
	}
	return Pos;
}
