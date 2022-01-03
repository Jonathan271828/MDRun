#ifndef _Forces
#include "force.h"
#endif
#ifndef _OUTPUT
#include "output.h"
#endif



void OnSitePotential::OnSitePotentialInit( const Real CP , const int Order ,
		const std::vector<VecNSlice<Real> >& PosEqui ){

	auto col = 0;
	for ( auto i = 0 ; i < PosEqui.size() ; i++ ){
		for ( auto j = 0 ; j < PosEqui[i].size() ; j++ ){
			EquiPosArray.push_back( PosEqui[i][j] );
			col++;
		}
	}

	unsigned int dim  = ( unsigned int )( PosEqui[0].size() );
	AssociateVecNSlice( EquiPosArray , EquiPos , dim );

	//Coupling.resize( Order - 1 );
	Coupling.push_back( CP );



	//for ( auto i = 0 ; i < Coupling.size() ; i++ ){
	//	if ( i%2 == 0 ){
	//		Coupling[i] = CP * Real( ( i+1 ) ) 
	//			       /  Real( Factorial( i+2 ));

	//	}
	//	else{
	//		Coupling[i] = -Real( ( i+1 ) ) /  Real( Factorial( i+2 ) );
	//	}
	//}
}



unsigned int OnSitePotential::Factorial( unsigned int N ){


	if ( N != 0 ){
		unsigned int Result = 1 ;
		for ( auto i = 1 ; i <= N ; i++ ){
			Result = Result * i;
		}
		return Result;
	}
	else{
		return 0;
	}
}



Real OnSitePotential::ComputeEnergyOnSite( const std::vector<VecNSlice<Real> >&
		                      Pos , 
				      const std::vector<std::vector<Real> >& Box ){

	Real TotE  =  Real( 0 );
#ifdef USE_OPENMP
#pragma omp parallel for reduction(+:TotE)
#endif
	for ( auto i = 0 ; i < Pos.size() ; i++ ){
		// get periodic image of EquiPos
		std::vector<Real> image = get_nearest_image( EquiPos[i] , Pos[i] );
		// compute vector point towards center
		image        =  EquiPos[i] - image;
		Real  norm   =  VecNorm( MatVec( Box , image ) );
		TotE   =  TotE  +  EnergySingleAtom( norm );
	}

	return TotE;
}



Real OnSitePotential::EnergySingleAtom( const Real norm ){


	Real Energy     =  Real( 0 );
	Real normPower  =  norm;
	for ( auto i = 0 ; i < Coupling.size() ; i++ ){
		normPower    =  normPower * normPower;
		Energy  =  Energy + Coupling[i] * normPower;
	}
	return Energy;
}


Real OnSitePotential::ComputeEnergySingleAtomOnSite( const VecNSlice<Real>& Pos ,
		                                     const std::vector<std::vector<Real> >&Box,
	     	                                     const int index ){


	std::vector<Real> image = get_nearest_image( EquiPos[index] , Pos );
	// compute vector point towards center
	image        =  EquiPos[index] - image;
	Real  norm   =  VecNorm( MatVec( Box , image ) );
	return EnergySingleAtom( norm );

}




void OnSitePotential::ComputeForceOnSite( const std::vector<VecNSlice<Real> > & Pos ,
		                          const std::vector<std::vector<Real> >& Box ,
                                          const std::vector<std::vector<Real> >& RecBox ,
		                          std::vector<VecNSlice<Real> >& Force ){
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
	for ( auto i = 0 ; i < Pos.size() ; i++ ){
		// get periodic image of EquiPos
		//EquiPos[i].Print();
		//Pos[i].Print();
		std::vector<Real> image = get_nearest_image( EquiPos[i] , Pos[i] );
		// compute vector point towards center
		image        =  MatVec( Box , EquiPos[i] - image );
		Real norm    =  VecNorm( image );
		Force[i]     =  0;
		if ( norm != Real( 0 ) ) image  =  image / norm;
		//std::cout << norm << std::endl;
	        for ( auto j = 0 ; j < Coupling.size() ; j++ ){
			norm  =  norm  *  norm;
			Force[i]  =  Force[i]  +  Coupling[j] * norm * image;
                }
		Force[i]  =  MatVec( RecBox , Force[i] );
        }
}


std::vector<Real> OnSitePotential::ComputeForceOnSiteSingle( 
		                       const VecNSlice<Real> & Pos ,
		                       const std::vector<std::vector<Real> >& Box ,
				       const std::vector<std::vector<Real> >& RecBox , 
				       const int index ){



	std::vector<Real> force( Pos.size() );

	std::vector<Real> image = get_nearest_image( EquiPos[index] , Pos );
	// compute vector point towards center
	image        =   MatVec( Box , EquiPos[ index ] - image );
	Real norm    =   VecNorm( image );
	if ( norm != Real( 0 ) ) image  =  image / norm;
	//std::cout << norm << std::endl;
	for ( auto j = 0 ; j < Coupling.size() ; j++ ){
		norm  =  norm  *  norm;
		force  =  force  +  Coupling[j] * norm * image;
        }
	force  =  MatVec( RecBox , force );
	return force;
}













/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 *
 */


 /*
  *
  *
  * double well coupling between atoms
  *
  *
  */

 void DoubleWellPotential::SetDoubleWellParameters( Real A , Real B , Real C , 
		                                    Real EquiDistCart ){
	 DoubleWellParA          =   A ;
	 DoubleWellParB          =   B ;
	 DoubleWellParCCart      =   C * A ;
	 printf( "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" );
	 printf( "Double Well potential initialized with \n" );
	 printf( "parameter A = %f and \n" , A );
	 printf( "parameter B = %f \n" , B );
	 printf( "parameter C = %f \n" , C * A );
	 printf( "Two particle double well potential initialized\n" );
	 printf( "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" );
	 DoubleWellEquiDistCart  =  EquiDistCart;
	 Normfactor = 1.0 / sqrt( M_PI / DoubleWellParCCart );
	 ExpFactor  = Real( 2.0 ) * Normfactor * DoubleWellParCCart * DoubleWellParB;
	 HarmFactor = Real( 2.0 ) * DoubleWellParA;
 }


std::vector<Real> DoubleWellPotential::DoubleWellForcePairPB( const VecNSlice<Real> & PosA , 
	  	                                              const VecNSlice<Real> & PosB ,
							      const std::vector<std::vector<Real>>& Box,
							      const std::vector<std::vector<Real>>& RecBox ){

	 std::vector<Real> r  =
	      get_nearest_image( PosA , PosB ) - PosA;
	 r           =   MatVec( Box , r );
	 Real norm   =   VecNorm( r );
	 Real Rnorm  =   Real( 1 ) / norm;
	 Real dist   =   ( norm - DoubleWellEquiDistCart );
         Real dist2  =   DoubleWellParCCart * dist * dist;

	 std::vector<Real> force = ( HarmFactor * dist - ExpFactor * dist * exp( -dist2 ) ) * Rnorm * r;
	 force  =  MatVec( RecBox , force );
	 return force;
 }


std::vector<Real> DoubleWellPotential::DoubleWellForcePair( const VecNSlice<Real> & PosA , 
		                                            const VecNSlice<Real> & PosB , 
							    const std::vector<std::vector<Real>>& Box,
							    const std::vector<std::vector<Real>>& RecBox ){

	 std::vector<Real> r  =  PosB  -  PosA;
	 r           =  MatVec( Box , r );
	 Real norm   =  VecNorm( r );
	 Real Rnorm  =  Real( 1 ) / norm;
	 Real dist   =  ( norm - DoubleWellEquiDistCart );
         Real dist2  =   DoubleWellParCCart * dist * dist;


	 std::vector<Real> force = ( HarmFactor * dist - ExpFactor * dist * exp( -dist2 ) ) * Rnorm * r;
	 force  =  MatVec( RecBox , force );
	 return force;
 }


Real DoubleWellPotential::DoubleWellEpotPairPB( const VecNSlice<Real>& PosA , const VecNSlice<Real>& PosB ,
	                                        const std::vector<std::vector<Real> >& Box ){
	 std::vector<Real> r =
		 get_nearest_image( PosA , PosB ) - PosA ;

	 Real norm     =   VecNorm( MatVec( Box , r ) );
	 Real dist     =   norm  -  DoubleWellEquiDistCart;
	 Real dist2    =   dist * dist;
	 Real dist2p   =   DoubleWellParCCart * dist2;
	 Real Epot     =   DoubleWellParA * dist2  +  Normfactor * DoubleWellParB * exp( -dist2p );
	 return Epot;
} 


Real DoubleWellPotential::DoubleWellEpotPair( const VecNSlice<Real>& PosA , const VecNSlice<Real>& PosB ,
		                              const std::vector<std::vector<Real> >& Box ){
	 std::vector<Real> r  =  PosB - PosA ;
	 Real norm     =   VecNorm( MatVec( Box , r ) );
	 Real dist     =   norm  -  DoubleWellEquiDistCart;
	 Real dist2    =   dist * dist;
	 Real dist2p   =   DoubleWellParCCart * dist2;
	 Real Epot     =   DoubleWellParA * dist2  +  Normfactor * DoubleWellParB * exp( -dist2p );
	 return Epot;
} 

 /*
  *
  *
  * harmonic coupling between atoms
  *
  *
  */

 void HarmonicForce::SetHarmonicForceParameters( Real x , Real EquiDistCar ){

	 HarmonicPotPar       =   x;
	 HarmonicEquiDistCar  =   EquiDistCar;

	 printf( "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" );
	 printf( "Harmonic potential initialized\n" );
	 printf( "Coupling constant set to %f eV/A \n" , HarmonicPotPar );
	 printf( "Equilibrium distance in cartesian coords %f\n" , HarmonicEquiDistCar );
	 printf( "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" );
 }



 // computes force between pair of atoms taking into account periodic boundray
 // conditions and allows to supply individual coupling
 // contant between considered pair
 // input are direct coords
 std::vector<Real> HarmonicForce::HarmonicForcePairPB( const VecNSlice<Real>& PosA ,
		                                       const VecNSlice<Real>& PosB , 
						       const std::vector<std::vector<Real>> & Box ,
						       const std::vector<std::vector<Real>> & RecBox ){

	std::vector<Real> r =
	      get_nearest_image( PosA , PosB );

	// cartesian distance vector
	r  =  MatVec( Box , r - PosA );
	
	Real norm =  VecNorm( r );
	Real Rnorm = Real( 1 ) / norm;

	// force in cartesian coords
	std::vector<Real> force =  Real( 2 ) * HarmonicPotPar *
		                      ( norm - HarmonicEquiDistCar ) * r * Rnorm;
	// force in direct coordinates
	force = MatVec( RecBox , force );	
	return force;
 }


 // routine assumes that posA and be are already given such that
 // a PosA and PosB are not minimum image convention and allows to supply
 // individual coupling contant between considered pair
 // supply positions in direct coordinates
 std::vector<Real> HarmonicForce::HarmonicForcePair( const VecNSlice<Real>& PosA , const VecNSlice<Real>& PosB ,
	        			             const std::vector<std::vector<Real>> & Box ,
						     const std::vector<std::vector<Real>> & RecBox ){
	 std::vector<Real> r = PosB - PosA ;
	 r = MatVec( Box , r );
	 Real norm  =  VecNorm( r );
	 Real Rnorm =  Real( 1 ) / norm;
	 std::vector<Real> force = Real( 2 ) * HarmonicPotPar * ( norm - HarmonicEquiDistCar ) * r * Rnorm;
	 force      =  MatVec( RecBox , force );
	 return force;
 }




 Real HarmonicForce::ComputePotential( const std::vector<VecNSlice<Real> >& Pos ,
                                       const std::vector<std::vector<int> >& LinkList ,
				       const std::vector<std::vector<Real> >& Box ){

	Real PotentialE = Real( 0.0 );
	for ( auto i = 0 ; i < Pos.size() ; i++ ){
		for ( auto j = 0 ; j < LinkList[i].size() ; j++ ){
			std::vector<Real> r =
				get_nearest_image( Pos[i] , Pos[LinkList[i][j]] ) - Pos[i];
		        Real norm = VecNorm( MatVec( Box , r ) );
			PotentialE = PotentialE  + HarmonicPotPar * ( norm - HarmonicEquiDistCar ) *
				                                    ( norm - HarmonicEquiDistCar );
	        }
        }
	return PotentialE;
 }

 /*
  * 2 potential energy functions are following that correspond exactly to the
  * above defined force functions
  */


 Real HarmonicForce::HarmonicEpotPair( std::vector<Real> PosA , std::vector<Real> PosB ,
		                       std::vector<std::vector<Real>> Box ){
	 std::vector<Real> r  =  PosB - PosA ;
	 Real norm = VecNorm( MatVec( Box , r ) );
	 Real Epot = HarmonicPotPar * ( norm - HarmonicEquiDistCar ) * ( norm - HarmonicEquiDistCar );
	 return Epot;
 }


 Real HarmonicForce::HarmonicEpotPairPB( const VecNSlice<Real>& PosA , const VecNSlice<Real>& PosB ,
		                                 const std::vector<std::vector<Real> >& Box ){
	 std::vector<Real> r =
		 get_nearest_image( PosA , PosB ) - PosA ;
	 Real norm = VecNorm( MatVec( Box , r ) );
	 Real Epot = HarmonicPotPar * ( norm - HarmonicEquiDistCar ) * ( norm - HarmonicEquiDistCar );
	 return Epot;
 }



/********************************************
 ********************************************
 ********************************************
 *
 * Morse potential coupling between atoms
 *
 ********************************************
 ********************************************
 ********************************************/


void MorsePotential::SetMorseForceParameters( Real x , Real descentCart , Real EquiDistCar ){

	MorseD            =   x;
	MorsealphaCart    =   descentCart;
	MorseEquiDistCart  =   EquiDistCar;
	printf( "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" );
	printf( "Morse potential initialized\n" );
	printf( "Potential depth set to %f eV \n" , MorseD );
	printf( "steepnes of potential set to %f [1/A]\n" , descentCart );
	printf( "Equilibrium distance in cartesian coords %f\n" , MorseEquiDistCart );
	printf( "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" );
}



// computes force between pair of atoms taking into account periodic boundray
 // conditions and allows to supply individual coupling
 // contant between considered pair
 std::vector<Real> MorsePotential::MorseForcePairPB( const VecNSlice<Real>& PosA , const VecNSlice<Real>& PosB ,
		                                     const std::vector<std::vector<Real> >& Box , 
		                                     const std::vector<std::vector<Real> >& RecBox ){

	 std::vector<Real> r =
	      get_nearest_image( PosA , PosB ) - PosA;
	 // cartesian connection vector
	 r  =  MatVec( Box , r );
	 // cartesian norm
	 Real norm = VecNorm( r );
	 Real Rnorm = Real( 1 ) / norm;
	 Real ExpFactor = exp( -MorsealphaCart * ( norm - MorseEquiDistCart ) );
	 ExpFactor = ( Real( 1.0 ) - ExpFactor ) * ExpFactor;
	 // force in cartesian coordinates
	 std::vector<Real> force = Real( 2.0 ) * MorseD * r * MorsealphaCart * Rnorm * ExpFactor;
	 // force in direct coordinates
	 force  =  MatVec( RecBox , force );
	 return force;
 }


 // routine assumes that posA and be are already given such that
 // a PosA and PosB are minimum image convention and allows to supply
 // individual coupling contant between considered pair
 std::vector<Real> MorsePotential::MorseForcePair( const VecNSlice<Real>& PosA , const VecNSlice<Real>& PosB, 
		                                   const std::vector<std::vector<Real> >& Box , 
		                                   const std::vector<std::vector<Real> >& RecBox ){


	 std::vector<Real> r = PosB - PosA ;
	 Real norm = VecNorm( r );
	 // cartesian connection vector
	 r  =  MatVec( Box , r );
	 // cartesian norm
	 Real Rnorm = Real( 1 ) / norm;
	 Real ExpFactor = exp( -MorsealphaCart * ( norm - MorseEquiDistCart ) );
	 ExpFactor = ( Real( 1.0 ) - ExpFactor ) * ExpFactor;
	 // force in cartesian coordinates
	 std::vector<Real> force = Real( 2.0 ) * MorseD * r * MorsealphaCart * Rnorm * ExpFactor;
	 // force in direct coordinates
	 force  =  MatVec( RecBox , force );
	 return force;
 }


 /*
  * 2 potential energy functions are following; that correspond exactly to the
  * above defined force functions with and without peridoc boundary condistions
  */


 Real MorsePotential::MorseEpotPair( const VecNSlice<Real>& PosA , const VecNSlice<Real>& PosB ,
	  	                     const std::vector<std::vector<Real> >& Box ){
	 std::vector<Real> r  =  PosB - PosA ;
	 Real norm = VecNorm( MatVec( Box , r ) );
	 Real ExpFactor = exp( -MorsealphaCart * ( norm - MorseEquiDistCart ) );
	 Real Epot = MorseD * ( Real( 1.0 ) - ExpFactor ) * ( Real( 1.0 ) - ExpFactor );
			   
	 return Epot;
 }



 Real MorsePotential::MorseEpotPairPB( const VecNSlice<Real>& PosA , const VecNSlice<Real>& PosB ,
 		                       const std::vector<std::vector<Real> >& Box ){
	 std::vector<Real> r =
		 get_nearest_image( PosA , PosB ) - PosA ;

	 Real norm = VecNorm( MatVec( Box , r ) );
	 Real ExpFactor = exp( -MorsealphaCart * ( norm - MorseEquiDistCart ) );
	 Real Epot = MorseD * ( Real( 1.0 ) - ExpFactor ) * ( Real( 1.0 ) - ExpFactor );
	 return Epot;
 }






 void ComputeForces::ComputeMorsePotentialEnergyArray( const std::vector<VecNSlice<Real> >& Pos , 
                                                        const std::vector<std::vector<Real> >& lattice ,
							std::vector<Real>& Epot ){

	 for ( size_t i = 0 ; i < Pos.size() ; i++ ){
		 Real energy = Real(0);
		 for ( size_t j = 0 ; j < LinkList.NNList[i].size() ; j++){
			 VecNSlice<Real> TempPos = Pos[ LinkList.NNList[i][j] ];
			 energy += MorseEpotPairPB( Pos[i] , TempPos , lattice ); 
		 }
		 Epot[i] = energy;
	 }
 }



 /*
  *
  *
  * class is inheriting from single force types
  * and computes total force on the considered atoms
  *
  *
  */

 ComputeForces::ComputeForces( const std::vector<VecNSlice<Real> >& Pos , unsigned int NN ,
		               const std::vector<std::vector<Real> >& Box ){
	 LinkList.ComputeNNListPeriodicBound( Pos , NN , Box );
 }


 // initialize force module with verlet cell list
 ComputeForces::ComputeForces( const std::vector<VecNSlice<Real> >& PosCar , const std::vector<VecNSlice<Real> >& PosDir ,
                               Real Cutoff , const std::vector<std::vector<Real> >& Box ){

	 LinkList.ComputeNNListPeriodicBoundInit( PosCar , PosDir , Cutoff , Box );
	 LinkList.PrintNNList();
 }



 void ComputeForces::SetUpOnSitePotential( const Real CP , const int Order ,  
		 const std::vector<VecNSlice<Real> >& PosEqui ){

	 OnSitePotentialInit( CP , Order , PosEqui );
 }



 //  initialize harmonic pair potential
 void ComputeForces::SetUpHarmonicPart( Real x , Real EquiDistCar ){
	 SetHarmonicForceParameters( x , EquiDistCar );
 }


 // initialize morse potential
 void ComputeForces::SetUpMorsePart( Real D , Real alphaCart , Real EquiDistCar ){
	 SetMorseForceParameters( D , alphaCart , EquiDistCar );
 }




 // initialize mexican hat onsite potential
 void ComputeForces::SetUpDoubleWellPart( Real A , Real B , Real C , 
		                          Real EquiDist ){
	 SetDoubleWellParameters( A , B , C , EquiDist );
 }


 //
 // only double well coupling between atoms
 void ComputeForces::SupplyDoubleWellOnly( const std::vector<VecNSlice<Real> >& Pos ,
		                           const std::vector<std::vector<Real> > & Box ,
		                           const std::vector<std::vector<Real> > & RecBox ,
   		                           std::vector<VecNSlice<Real> >& force ){

#ifdef USE_OPENMP
#pragma omp parallel for
#endif
	for ( auto i = 0 ; i < Pos.size() ; i++ ){
		force[i] = 0;
		for ( auto j = 0 ; j < LinkList.NNList[i].size() ; j++ ){
			force[i] = force[i] + DoubleWellForcePairPB( Pos[i] , Pos[ LinkList.NNList[i][j]] ,
			  		                             Box , RecBox );
		}
	}
 }


 Real ComputeForces::PotentialEnergyDoubleWell( const std::vector<VecNSlice<Real> >& Pos,
		                                const std::vector<std::vector<Real> >& Box ){
	 Real Epot = Real( 0 );
#ifdef USE_OPENMP
#pragma omp parallel for reduction(+:Epot)
#endif
	 for ( auto i = 0 ; i < Pos.size() ; i++ ){
		 for ( auto j = 0 ; j < LinkList.NNList[i].size() ; j++ ){
		         Epot = Epot  +  DoubleWellEpotPairPB( Pos[i] , Pos[ LinkList.NNList[i][j]] , Box );
		 }
         }
         return Epot;
 }
 // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~




 //
 // only harmonic coupling between atoms
 // input are direct coordinates
 // output is force in direct coordinates
 void ComputeForces::SupplyHarmonicOnly( const std::vector<VecNSlice<Real> >& Pos ,
                                         const std::vector<std::vector<Real> >& Box ,
                                         const std::vector<std::vector<Real> >& RecBox ,
   		                         std::vector<VecNSlice<Real> >& force ){

#ifdef USE_OPENMP
#pragma omp parallel for
#endif
	for ( auto i = 0 ; i < Pos.size() ; i++ ){
		force[i] = 0;
		for ( auto j = 0 ; j < LinkList.NNList[i].size() ; j++ ){
			force[i] = force[i] + HarmonicForcePairPB( Pos[i] , Pos[ LinkList.NNList[i][j]] ,
					                Box , RecBox );
		}
	}
 }


 Real ComputeForces::PotentialEnergyHarmonic( const std::vector<VecNSlice<Real> >& Pos,
		                              const std::vector<std::vector<Real> >& Box ){
	 Real Epot = Real( 0 );
#ifdef USE_OPENMP
#pragma omp parallel for reduction(+:Epot)
#endif
	 for ( auto i = 0 ; i < Pos.size() ; i++ ){
		 for ( auto j = 0 ; j < LinkList.NNList[i].size() ; j++ ){
		         Epot = Epot + HarmonicEpotPairPB( Pos[i] , Pos[ LinkList.NNList[i][j]] , Box );
		 }
         }
         return Epot;
 }
 // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

 //
 // only Morse potential like coupling between atoms
 //
 // function computes morse potential for certain atom pair-> note those atom pair has to be nearest neighbour
 // otherwise it does not make sense can be called from outside too
 Real ComputeForces::SupplyMorsePairOnlyPB( const VecNSlice<Real>& PosA , const VecNSlice<Real>& PosB , 
		                                    const std::vector<std::vector<Real > >& Box ){
	 return MorseEpotPairPB( PosA , PosB , Box );
 }

 // computes the energy for atom with index AtNr
 Real ComputeForces::SupplyMorseNNEOnlyPB( unsigned int AtNr , const std::vector<VecNSlice<Real> >& Pos , 
		                                   const std::vector<std::vector<Real> >& Box ){
	 Real Epot = Real( 0.0 );
	 for ( auto i = 0 ; i < LinkList.NNList[AtNr].size() ; i++ ){
		 Epot = Epot + MorseEpotPairPB( Pos[AtNr] , Pos[ LinkList.NNList[AtNr][i]] , Box );
	 }
	 return Epot;
 }

 // compute potential energy for certain atom with neighbor list List
 Real ComputeForces::SupplyMorseSingleAtom( const VecNSlice<Real>& PosNew , 
		                            const std::vector<VecNSlice<Real> >& Pos , 
		                            const std::vector<unsigned int>& List ,
					    const std::vector<std::vector<Real> >& Box ){
	 Real Epot = Real( 0.0 );
	 for ( auto i = 0 ; i < List.size() ; i++ ){
		 Epot = Epot + MorseEpotPairPB( PosNew , Pos[ List[i] ] , Box );
	 }
	 return Epot;
 }


 // compute potential force for certain atom with neighbor list List
std::vector<Real> ComputeForces::SupplyMorseForceSingleAtom( const VecNSlice<Real>& PosNew ,
		                                             const std::vector<VecNSlice<Real> >& Pos , 
		                                             const std::vector<unsigned int>& List ,
							     const std::vector<std::vector<Real>>& Box ,
							     const std::vector<std::vector<Real>>& RecBox ){

	 std::vector<Real> force ;
	 force.resize( PosNew.size() );
	 for ( auto i = 0 ; i < List.size() ; i++ ){
		 force = force + MorseForcePairPB( PosNew , Pos[ List[i] ] , Box , RecBox );
	 }
	 return force;
 }




 void ComputeForces::SupplyMorseOnly( const std::vector<VecNSlice<Real> >& Pos , 
		                      const std::vector<std::vector<Real> > & Box ,
		                      const std::vector<std::vector<Real> > & RecBox ,
		                      std::vector<VecNSlice<Real> >& force ){

#ifdef USE_OPENMP
#pragma omp parallel for
#endif
	for ( auto i = 0 ; i < Pos.size() ; i++ ){
		//rezero force components
		for ( auto j = 0 ; j < Pos[i].size() ; j ++ ){
			force[i][j]  =  Real( 0 );
		}
		// compute pair interactions and add them up
		for ( auto j = 0 ; j < LinkList.NNList[i].size() ; j ++ ){
			force[i] = force[i] + MorseForcePairPB( Pos[i] , Pos[ LinkList.NNList[i][j]] ,
					                        Box , RecBox );
		}
	}
}

 Real ComputeForces::PotentialEnergyMorse( const std::vector<VecNSlice<Real> >& Pos , 
		                           const std::vector<std::vector<Real> >& Box ){

	 Real Epot = Real( 0 );
#ifdef USE_OPENMP
#pragma omp parallel for reduction(+:Epot)
#endif
	 for ( auto i = 0 ; i < Pos.size() ; i++ ){
		 for ( auto j = 0 ; j < LinkList.NNList[i].size() ; j++ ){
		         Epot = Epot + MorseEpotPairPB( Pos[i] , Pos[ LinkList.NNList[i][j] ] , Box );
		 }
     }
	 return Epot;
 }

 // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



 //
 // compute harmonic coupling between atoms and constrain atoms to their crystal site by harmonic potential
// void ComputeForces::SupplyHarmonicAndOnSite( const std::vector<VecNSlice<Real> >& Pos ,
//		                                    std::vector<VecNSlice<Real> >& force ){
//
//	for ( auto i = 0 ; i < Pos.size() ; i++ ){
//		for ( auto j = 0 ; j < LinkList.NNList[i].size() ; j ++ ){
//			force[i] = force[i] + HarmonicForcePairPB( Pos[i] , Pos[ LinkList.NNList[i][j]] );
//		}
//	}
// }

 Real ComputeForces::PotentialEnergyHarmonicAndDoubleWell( const std::vector<VecNSlice<Real> >& Pos , 
		                                           const std::vector<std::vector<Real> >& Box ){


	 Real Epot = Real( 0 );
#ifdef USE_OPENMP
#pragma omp parallel for reduction(+:Epot) 
#endif
	 for ( auto i = 0 ; i < Pos.size() ; i++ ){
		 for ( auto j = 0 ; j < LinkList.NNList[i].size() ; j++ ){
		         Epot = Epot  +  HarmonicEpotPairPB( Pos[i] , Pos[ LinkList.NNList[i][j]] , Box );
		         Epot = Epot  +  DoubleWellEpotPairPB( Pos[i] , Pos[ LinkList.NNList[i][j]] , Box );
		 }
     }
     return Epot;
 }



 void ComputeForces::SupplyOnSiteForceOnly( const std::vector<VecNSlice<Real> >& Pos ,
		                            const std::vector<std::vector<Real> >& Box ,
					    const std::vector<std::vector<Real> > &RecBox ,
		                            std::vector<VecNSlice<Real> >& force ){
	 
	 ComputeForceOnSite( Pos , Box , RecBox , force );
 }


 Real ComputeForces::SupplyOnSiteEnergyOnly( const std::vector<VecNSlice<Real> >& Pos ,
		                             const std::vector<std::vector<Real> >& Box ){
	 return ComputeEnergyOnSite( Pos , Box );
 }




 void ComputeForces::SupplyHarmonicAndOnSiteForce( const std::vector<VecNSlice<Real> >& Pos ,
		                                   const std::vector<std::vector<Real> >& Box ,
					           const std::vector<std::vector<Real> > &RecBox ,
						   std::vector<VecNSlice<Real> >& force ){

#ifdef USE_OPENMP
#pragma omp parallel for 
#endif
	 for ( auto i = 0 ; i < Pos.size() ; i++ ){
		 force[i]  =  ComputeForceOnSiteSingle( Pos[i] , Box , RecBox , i );
		 for ( auto j = 0 ; j < LinkList.NNList[i].size() ; j++ ){
			 force[i]  =  force[i]  + HarmonicForcePair( Pos[i] , 
					             Pos[LinkList.NNList[i][j]] , Box , RecBox ); 
		 }
	 }
 }


 Real ComputeForces::PotentialEnergyHarmonicAndOnSite( const std::vector<VecNSlice<Real> >& Pos ,
	           	              const std::vector<std::vector<Real> >& Box ){

	 Real Epot = Real( 0 );
#ifdef USE_OPENMP
#pragma omp parallel for reduction(+:Epot) 
#endif
	 for ( auto i = 0 ; i < Pos.size() ; i++ ){
		 Epot  =  Epot  =  ComputeEnergySingleAtomOnSite( Pos[i] , Box , i );
		 for ( auto j = 0 ; j < LinkList.NNList[i].size() ; j++ ){
		         Epot = Epot  +  HarmonicEpotPairPB( Pos[i] , Pos[ LinkList.NNList[i][j]] , Box );
		 }
	 }

	 return Epot;
 }


