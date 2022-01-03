
#include "thermostats.h"






 /*
  * initialize the thermostat with temperature
  * collison probability nu and the sigma
  * for the distribution from which the
  * velocities are taken
  *
  */

 void ThermostatsVeloVerlet::InitAndersen( const Real temp , const Real nu , const std::vector<Real>& mass ){


	 Sigma.resize( mass.size() );
	 for ( auto i = 0 ; i < mass.size(); i++ ){
		 Sigma[i] = sqrt( units.kB * temp /mass[i] );
	 }
	 ColProb = nu;
	 ThermoType = 0;
 }


 // Nose-Hoover initialization
 void ThermostatsVeloVerlet::InitNoseHoover( Real q , unsigned int DOF , Real temp , const std::vector<Real>& mass ){

	 VirtMass = Real( DOF ) * temp * units.kB * q*q ;
	 ThermoType = 1;
	 NoseG = Real( DOF );
	 temperature = temp*units.kB;
	 Sigma.resize( mass.size() );
	 for ( auto i = 0 ; i < mass.size(); i++ ){
		 Sigma[i] = sqrt( units.kB * temp /mass[i] );
	 }
 }


 // Langevin initialization
 void ThermostatsVeloVerlet::InitLangevin( Real Fric , Real temp , Real dt ,
		                           unsigned int dim , const std::vector<Real>& Mass ){

	 LangevinFric = Fric ;
	 ThermoType = 2;
	 temperature = temp*units.kB;

	 LangevinA = ( Real( 2 ) - Fric*dt ) / ( Real( 2 ) + Fric*dt );
	 LangevinC = Real( 2 ) * dt / ( Real( 2 ) + Fric * dt );


	 Sigma.resize( Mass.size() );
	 RandomForceData.resize( Mass.size() * dim );
	 LangevinB.resize( Mass.size() );
	 for ( auto i = 0 ; i < Mass.size(); i++ ){
		 Sigma[i] = sqrt( units.kB * temp / Mass[i] );
	     LangevinB[i] = sqrt( temperature*Fric * Real( 0.5 ) * dt / Mass[i] );
	 }


	 for ( auto i = 0 ; i < RandomForceData.size() ; i=i+dim ){
		 VecNSlice<Real> A( RandomForceData.data() + i , dim );
		 RandomForce.push_back( A );
	 }
	 LangevinRandomNumbers();
 }

// ~?~?~?~?~?~?~?~?~?~?~?~?~?~?~?~?~?~?~?~?~?~?~?~?~?~?~?~?

 /// routine to initialize the Bussi Donadio Parrinello thermostat
 // implemented according to Canonical sampling through velocity
 // rescaling
 // The Journal of Chemical Physics 126, 014101 (2007)
void ThermostatsVeloVerlet::InitBDP( Real tau , Real Temp , Real dt , int DOF , 
   		                     const std::vector<Real>& mass ){

	BDPKinetic  =   Real( DOF ) * units.kB * Temp / Real( 2.0 );
	BDPFactor1  =   exp( -dt / tau );
	BDPFactor2  =   2 * exp( -0.5 * dt / tau );
	BDPFactor3  =   1.0 - BDPFactor1;
	BDPDOF      =   DOF;
	if ( ( DOF - 1 )% 2 == 0 ){
		evenDOF = false;
	        BDPDOFMOD2  =   ( DOF - 1 ) ;
	}
	else{
		// take one less DOF 
		evenDOF = true;
	        BDPDOFMOD2  =   ( DOF - 2 ) ;
	}
	
	// to initialize velocisties
	Sigma.resize( mass.size() );
	for ( auto i = 0 ; i < mass.size(); i++ ){
	        Sigma[i] = sqrt( units.kB * Temp / mass[i] );
	}
	ThermoType = 3;
	BDPHeatGrad = false;
}


// helper routines for Bussi Donadio Parrinello thermostat
void ThermostatsVeloVerlet::BDPComputeKineticEnergy( const std::vector<VecNSlice<Real> >& VeloCity ,
		                        const std::vector<std::vector<Real> >&Box , 
					const std::vector<Real> & Mass ){

	Real Ek = Real( 0.0 );
#if USE_OPENMP
#pragma omp parallel for reduction(+:Ek)
#endif
	for ( auto i = 0; i < VeloCity.size() ; i++ ){
		std::vector<Real> Cart = MatVec( Box , VeloCity[i] );
		Ek = Ek  +  DotProduct( Cart , Cart ) * Mass[i];
	}
	BDPKineticAct = Ek * 0.5 ;
}

// helper routines for Bussi Donadio Parrinello thermostat
void ThermostatsVeloVerlet::BDPComputeKineticEnergyHeatGrad( 
		                        const std::vector<VecNSlice<Real> >& VeloCity ,
		                        const std::vector<std::vector<Real> >&Box , 
					const std::vector<Real> & Mass ){

#if USE_OPENMP
#pragma omp parallel for
#endif
	for ( auto i = 0 ; i < HeatGradLink.size(); i++ ){
		BDPKineticLayerAct[i] = Real( 0 );
		for ( auto j = 0 ; j < HeatGradLink[i].size() ;j++ ){ 
			auto indx = HeatGradLink[i][j];
			std::vector<Real> Cart = MatVec( Box , VeloCity[indx] );
			BDPKineticLayerAct[i] += DotProduct( Cart , Cart ) * Mass[indx];
		}
		BDPKineticLayerAct[i] = 0.5 * BDPKineticLayerAct[i];
	}
}




void ThermostatsVeloVerlet::BDPRecomputeAlpha( void ){


	Real Gr = Rand.generate_gamma_dist( BDPDOFMOD2 ); // gamma random var
	Real r  = Rand.generate_gauss_dist( 0.0 , 1.0 ); // gaussian radom var
	Real r2 = r*r;
	// adapt gamma distributed according to
	// uneven number of degrees of freedom
	if ( evenDOF ){
		Real x = Rand.generate_gauss_dist( 0.0 , 1.0 );
		Gr = Gr + x*x;
	}

	Gr = Gr + r2;


	Real factor = BDPKinetic / ( BDPDOF * BDPKineticAct ); // kinetic energy ratio

	BDPAlpha = BDPFactor1 + factor * BDPFactor3 * Gr 
		              + BDPFactor2 * sqrt( factor * BDPFactor3 ) * r;
	BDPAlpha = sqrt( BDPAlpha );
}


void ThermostatsVeloVerlet::BDPRecomputeAlphaHeatGrad( void ){


	for ( auto i = 0 ; i < HeatGradLink.size(); i++ ){
		Real Gr = Rand.generate_gamma_dist( BDPDOFMOD2 ); // gamma random var
	        Real r  = Rand.generate_gauss_dist( 0.0 , 1.0 ); // gaussian radom var
	        Real r2 = r*r;
	        // adapt gamma distributed according to
	        // uneven number of degrees of freedom
	        if ( evenDOF ){
	        	Real x = Rand.generate_gauss_dist( 0.0 , 1.0 );
	        	Gr = Gr + x*x;
	        }

	        Gr = Gr + r2;


		// kinetic energy ration in layer
	        Real factor = BDPKineticLayer[i] / ( BDPDOF * BDPKineticLayerAct[i] ); // kinetic energy ratio

	        BDPLayerAlpha[i] = BDPFactor1 + factor * BDPFactor3 * Gr 
	        	              + BDPFactor2 * sqrt( factor * BDPFactor3 ) * r;
	        BDPLayerAlpha[i] = sqrt( BDPLayerAlpha[i] );
        }
}


void ThermostatsVeloVerlet::BDPRecomputeAlphaHeatGradAM( void ){


	for ( auto i = 0 ; i < HeatGradLink.size(); i++ ){
		Real Gr = Rand.generate_gamma_dist( BDPDOFMOD2 ); // gamma random var
	        Real r  = Rand.generate_gauss_dist( 0.0 , 1.0 ); // gaussian radom var
	        Real r2 = r*r;
	        // adapt gamma distributed according to
	        // uneven number of degrees of freedom
	        if ( evenDOF ){
	        	Real x = Rand.generate_gauss_dist( 0.0 , 1.0 );
	        	Gr = Gr + x*x;
	        }

	        Gr = Gr + r2;


		// kinetic energy ration in layer
	        Real factor = BDPKineticLayer[i] / ( BDPDOF * BDPKineticLayerAct[i] ); // kinetic energy ratio

	        BDPLayerAlpha[i] = BDPFactor1 + factor * BDPFactor3 * Gr 
	        	              + BDPFactor2 * sqrt( factor * BDPFactor3 ) * r;
	        BDPLayerAlpha[i] = sqrt( BDPLayerAlpha[i] );
	}
}



// rescales velocities according to the computed alpha
void ThermostatsVeloVerlet::BDPRescaleVelocities( std::vector<VecNSlice< Real> >& VeloCity ){

	for ( auto i = 0 ; i < VeloCity.size(); i++ ){
		VeloCity[i]  =  VeloCity[i] * BDPAlpha;
	}
}

// rescales velocities for different layers of temperature
// only rescale the velocity in certain areas to get a real heat gradient
void ThermostatsVeloVerlet::BDPRescaleVelocitiesHeatGradAM( std::vector<VecNSlice< Real> >& VeloCity ){

	for ( auto i = 0 ; i < AMLinkList.size(); i++ ){
		auto indA = AMLinkList[ i ];
		for ( auto j = 0 ; j < HeatGradLink[indA].size(); j++ ){
			auto indx = HeatGradLink[ indA ][ j ];
	                VeloCity[ indx ] = VeloCity[ indx ] * BDPLayerAlpha[ indA ];
	        }
	}
}


// rescales velocities for different layers of temperature
void ThermostatsVeloVerlet::BDPRescaleVelocitiesHeatGrad( std::vector<VecNSlice< Real> >& VeloCity ){

	for ( auto i = 0 ; i < HeatGradLink.size(); i++ ){
		for ( auto j = 0 ; j < HeatGradLink[i].size(); j++ ){
			auto indx = HeatGradLink[ i ][ j ];
	                VeloCity[ indx ] = VeloCity[ indx ] * BDPLayerAlpha[ i ];
	        }
	}
}





// end helper routines Bussi Donadio Parrinello thermostat
// ~?~?~?~?~?~?~?~?~?~?~?~?~?~?~?~?~?~?~?~?~?~?~?~?~?~?~?~?

 
// heat gradient thermostating
// andersen heat gradient initialization
void ThermostatsVeloVerlet::InitHeatGradientThermoAndersen( const HeatGradientStruc& InData  ,
		                            const std::vector<int>& NDirs ,
		                            const std::vector<Real>& mass ){




	// computing the heat gradient
	HeatVars =  InData;
	InitHeatGradientLinkList( NDirs );
	InitHeatGradientTemperatures( NDirs );
	Sigma.resize( mass.size() );
	for ( auto i = 0 ; i < HeatGradLink.size(); i++ ){
		for ( auto j = 0 ; j < HeatGradLink[i].size(); j++ ){
			auto indx = HeatGradLink[i][j];
	                Sigma[ indx ] = sqrt( units.kB * TempGradient[i] /mass[indx] );
	        }
	}
	ColProb = HeatVars.nu;
	ThermoType = 0;

	// make temperature output
//	DebugFile = fopen( "TempGrad.out" , "w" );
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


 void ThermostatsVeloVerlet::InitHeatGradientThermoLangevin( const HeatGradientStruc& InData ,
		       Real Fric , Real dt , unsigned int dim , const std::vector<int>& NDirs , 
		          const std::vector<Real>& Mass ){

	 HeatVars =  InData;
	 LangevinFric = Fric ;
	 ThermoType = 2;

	 LangevinA = ( Real( 2 ) - Fric*dt ) / ( Real( 2 ) + Fric*dt );
	 LangevinC = Real( 2 ) * dt / ( Real( 2 ) + Fric * dt );



	 InitHeatGradientLinkList( NDirs );
	 InitHeatGradientTemperatures( NDirs );

	 Sigma.resize( Mass.size() );
	 RandomForceData.resize( Mass.size() * dim );
	 LangevinB.resize( Mass.size() );

	 for ( auto i = 0 ; i < HeatGradLink.size(); i++ ){
	 	for ( auto j = 0 ; j < HeatGradLink[i].size(); j++ ){
	 		auto indx = HeatGradLink[i][j];
	                 Sigma[ indx ] = sqrt( units.kB * TempGradient[i] / Mass[indx] );
	                 LangevinB[indx] = sqrt( units.kB * TempGradient[i] *Fric * Real( 0.5 ) * dt / Mass[indx] );
	         }
	 }

	 for ( auto i = 0 ; i < RandomForceData.size() ; i=i+dim ){
		 VecNSlice<Real> A( RandomForceData.data() + i , dim );
		 RandomForce.push_back( A );
	 }
	 LangevinRandomNumbers();
	 //DebugFile = fopen( "TempGrad.out" , "w" );
 }

 void ThermostatsVeloVerlet::InitHeatGradientThermoBDP( const HeatGradientStruc& InData , Real tau , Real dt , 
		                                        const std::vector<Real>& mass , const std::vector<int>& NDirs , 
							int dim ){



         HeatVars = InData;
	 InitHeatGradientLinkList( NDirs );
	 InitHeatGradientTemperatures( NDirs );
	 
	 int DOF = HeatGradLink[0].size() * dim;  // degrees of freedom per layer

	 Sigma.resize( mass.size() );
	 BDPKineticLayer.resize( HeatGradLink.size() );
	 BDPKineticLayerAct.resize( HeatGradLink.size() );
	 BDPLayerAlpha.resize( HeatGradLink.size() );
	 // to properly initialize velocities
	 for ( auto i = 0 ; i < HeatGradLink.size(); i++ ){
	  	BDPKineticLayer[ i ] = Real( DOF ) * units.kB * TempGradient[i] / Real( 2.0 );
	 	for ( auto j = 0 ; j < HeatGradLink[i].size(); j++ ){
	 		auto indx = HeatGradLink[i][j];
	                 Sigma[ indx ] = sqrt( units.kB * TempGradient[i] / mass[indx] );
	        }
	 }

 	 BDPFactor1  =   exp( -dt / tau );
 	 BDPFactor2  =   2 * exp( -0.5 * dt / tau );
 	 BDPFactor3  =   1.0 - BDPFactor1;
 	 BDPDOF      =   DOF;
 	 if ( ( DOF - 1 )% 2 == 0 ){
 	 	 evenDOF = false;
 	         BDPDOFMOD2  =   ( DOF - 1 ) ;
 	 }
 	 else{
 	 	 // take one less DOF 
 	  	 evenDOF = true;
 	         BDPDOFMOD2  =   ( DOF - 2 ) ;
 	 } 
 	 ThermoType = 3;
	 BDPHeatGrad = true;
	 DebugFile = fopen( "TempGrad.out" , "w" );
         
	 // always take edge layers


	 if ( HeatVars.ActiveMeasure == 1 ){
	      AMLinkList.push_back( 0 );
	      if ( HeatGradLink.size() % 2 == 0 ){
	              AMLinkList.push_back( HeatGradLink.size()/2-1 );
	              AMLinkList.push_back( HeatGradLink.size()/2 );
	      }
	      else{
	              AMLinkList.push_back( (AMLinkList.size()-1) / 2 );
	      }
	      AMLinkList.push_back( HeatGradLink.size()-1 );
	 }
	 EquiSteps = 0;
}


/// heat gradient helper routines
//
// linking atoms to temperature layer
void ThermostatsVeloVerlet::InitHeatGradientLinkList( const std::vector<int>& NDirs ){
	

	size_t col = 0;
	std::vector<int> Use;
	Use.resize( NDirs.size() );
	int HDir = HeatVars.Direction - ( int ) ( 1 );
	HeatGradLink.resize( NDirs[ HDir ] );
	for ( auto i = 0 ; i < NDirs[0] ; i++ ){
		Use[0] = i;
		for ( auto j = 0 ; j < NDirs[1] ; j++ ){
		        Use[1] = j;
			for ( auto k = 0 ; k < NDirs[2] ; k++ ){
		                Use[2] = k;
				HeatGradLink[ Use[ HDir ] ].push_back( col );
				col++;
			}
		}
	}
}


// compute target temperaturs in different layers
void ThermostatsVeloVerlet::InitHeatGradientTemperatures( const std::vector<int>& NDirs ){

	int HDir = HeatVars.Direction - ( int ) ( 1 );

	Real dT  =  Real( 2 ) * ( HeatVars.temperatureA - HeatVars.temperatureB ) / 
		    Real( NDirs[ HDir ] );
	TempGradient.resize( NDirs[ HDir ] );
	int N   =  floor( NDirs[ HDir ] / int( 2 ) );
	int Max =  NDirs[ HDir ];
	
	if ( Max % 2 == 0 ){
	    for ( auto i = 0 ; i < N ; i++ ){
	    	Real T = HeatVars.temperatureB + Real( i ) * dT + dT / 2.0;
	    	TempGradient[ N - i - 1 ] =  T ;
	    	TempGradient[ N + i ] =  T ;
	    }
	}
	else{
	    TempGradient[ N ] =  HeatVars.temperatureB + dT / 2.0;
	    for ( auto i = 0 ; i < N ; i++ ){
	    	Real T = HeatVars.temperatureB + Real( i ) * dT + dT / 2.0;
	    	TempGradient[ N - i - 1 ] =  T ;
	    	TempGradient[ N + i + 1 ] =  T ;
	    
	    }
	}
}


//
// Write Heat Gradient Link List to screen
//
void ThermostatsVeloVerlet::WriteHeatGradientLinkList( void ){
	for ( auto i = 0 ; i < HeatGradLink.size() ; i++ ){
		for ( auto j = 0 ; j < HeatGradLink[i].size() ; j++ ){
			std::cout << HeatGradLink[ i ][ j ] << "   ";
		}
		std::cout << std::endl;
		std::cout << std::endl;
	}
}


void ThermostatsVeloVerlet::ComputeLayerTemperature( const std::vector<VecNSlice<Real> >& Velos , 
			      const std::vector<Real>& Mass ){

	//fprintf( DebugFile , " %i  " , 0 );
	for ( auto i = 0 ; i < HeatGradLink.size() ; i++ ){
		Real Ek = 0;
		auto Nf = HeatGradLink[i].size();
		for ( auto j = 0 ; j < HeatGradLink[i].size() ; j++ ){
			auto indx = HeatGradLink[i][j];
			Ek += DotProduct( Velos[indx] , Velos[indx] ) * Mass[indx];
		}
		Ek *= 0.5;
		Real T = Real( 2 ) * Ek / Real( 3 ) / Real( Nf ) / units.kB;
		fprintf( DebugFile , " %f  " , T );
	}
	fprintf( DebugFile , " %s" , "\n" );
}

// end of heat gradient helper routines
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

void ThermostatsVeloVerlet::VelocityExchange( std::vector<VecNSlice<Real> >& Velos ){

	int N = 1;

	std::vector<Real> High( N , 0 );
	std::vector<int> HighIndx( N , 0 );
	std::vector<Real> Low( N , Real(1000000) );
	std::vector<int> LowIndx( N , 0 );
	
	for ( auto i = 0 ; i < Velos.size() ; i++ ){
		Real norm = VecNSliceNorm( Velos[i] );
		for ( auto j = 0 ; j < N ; j++ ){
			if ( norm > High[j] ){
				High[j]     =  norm;
				HighIndx[j] =  i;
				break;
			}
		}
		for ( auto j = 0 ; j < N ; j++ ){
			if ( norm < Low[j] ){
			        Low[j]      =  norm;
				LowIndx[j] =  i;
				break;
			}
		}
	}

	//for ( auto i = 0 ; i < N ; i++ ){
	//	std::cout << High[ i ] << "  " 
	//		  << Low[ i ] << std::endl;
	//}



	for ( auto i = 0 ; i < N ; i++ ){
		int Hindx = HighIndx[i];
		int Lindx = LowIndx[i];
		Real HNorm = VecNSliceNorm( Velos[Hindx] );
		Real LNorm = VecNSliceNorm( Velos[Lindx] );
		//std::cout << Low[i] / HNorm << "  " << High[i] / LNorm << std::endl;
		Velos[Hindx] = Velos[Hindx] / HNorm * Low[i];
		Velos[Lindx] = Velos[Lindx] / LNorm * High[i];
	}
}


// velocity verlet in kick drift form implementation
void ThermostatsVeloVerlet::LeapFrogStepA( std::vector<VecNSlice<Real> >& Pos , 
		                           std::vector<VecNSlice<Real> >& Velo , 
				           const std::vector<VecNSlice<Real> >& Force ,
				           const std::vector<Real>& Mass , Real dt ){
#if USE_OPENMP
#pragma omp parallel for
#endif
	for ( auto i = 0 ; i < Pos.size(); i++ ){
		std::vector<Real> df = Real( 0.5 ) * dt * ( Force[i] / Mass[i] );
		Velo[i] = Velo[i] + df;
		std::vector<Real> dv =  dt * Velo[i] ;
		Pos[i] = Pos[i] + dv;
		ApplyPeriodicBound( Pos[i] );
        }
}


void ThermostatsVeloVerlet::LeapFrogStepB( std::vector<VecNSlice<Real> >& Velo ,
		                           const std::vector<VecNSlice<Real> >& Force ,
		                           const std::vector<Real>& Mass ,
					   Real dt ){
#if USE_OPENMP
#pragma omp parallel for
#endif
	for ( auto i = 0 ; i < Velo.size() ; i++ ){
		std::vector<Real> df = Real( 0.5 ) * dt * ( Force[i] / Mass[i] );
		Velo[i] = Velo[i] + df;
	}
}





/// andersen implementation
void ThermostatsVeloVerlet::AndersenStepA( std::vector<VecNSlice<Real> >& Pos ,
		                           const std::vector<VecNSlice<Real> >& Force ,
					   std::vector<VecNSlice<Real> >& Velo , Real dt , 
					   const std::vector<Real>& Mass ){

#if USE_OPENMP
#pragma omp parallel for
#endif
	for ( auto i = 0 ; i < Pos.size(); i++ ){
		std::vector<Real> df = Real( 0.5 ) * dt * ( Force[i] / Mass[i] );
		Velo[i] = Velo[i] + df;
		std::vector<Real> dv =  dt * Velo[i] ;
		Pos[i] = Pos[i] + dv;
		ApplyPeriodicBound( Pos[i] );
        }
}


void ThermostatsVeloVerlet::AndersenStepB( std::vector<VecNSlice<Real> >& Velo ,
		                           const std::vector<VecNSlice<Real> >& Force ,
					   Real dt , const std::vector<std::vector<Real> >& Box ,
		                           const std::vector<Real>& Mass ){

#if USE_OPENMP
#pragma omp parallel for
#endif
	for ( auto i = 0 ; i < Velo.size() ; i++ ){
		Real r = Rand.doub();
		if ( r < ColProb ){
			for ( auto j = 0 ; j < Velo[i].size(); j++ ){
				Velo[i][j] = Rand.generate_gauss_dist( Real( 0 ) , Sigma[i] );
		        }
			Velo[i] = MatVec( Box , Velo[i] );
		}
		else{
		    std::vector<Real> df = Real( 0.5 ) * dt * ( Force[i] / Mass[i] );
		    Velo[i] = Velo[i] + df;
		}
	}
}
/// end andersen implementation
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 void ThermostatsVeloVerlet::InitVelocities( std::vector<VecNSlice<Real> >& Velos , int dim ,
		const std::vector<std::vector<Real> >& Box , const std::vector<Real>& Mass ){

	 Real Ek = 0 ;
	 for ( auto i = 0 ; i < Velos.size() ; i++ ){
		for ( auto j = 0 ; j < Velos[i].size() ; j++ ){
			 Velos[i][j]  =  Rand.generate_gauss_dist( Real( 0 ) , Sigma[i] );
		}
		if ( ThermoType == 1 ){
			std::vector<Real> v = MatVec( Box , Velos[i] );
			Ek = Ek + DotProduct( v , v ) * Mass[i];
		}
	 }
	 if ( ThermoType == 1 ){
		 NoseChi = ( Ek - NoseG * temperature ) / VirtMass;
	 }
 }




 /// Nose Hoover implementation
 void ThermostatsVeloVerlet::NoseHooverStepA( std::vector<VecNSlice<Real> >& Pos , 
		                        const std::vector<VecNSlice<Real> >& Force , 
	      		                std::vector<VecNSlice<Real> >& Velo  , 
					Real dt , 
					const std::vector<std::vector<Real> >& Box , 
					const std::vector<Real>& Mass ){
	 Real Ek = Real( 0 );
#if USE_OPENMP
#pragma omp parallel for reduction(+:Ek)
#endif
	 for ( auto i = 0 ; i < Pos.size() ; i++ ){
		 std::vector<Real> dv = dt * Velo[i];
		 //std::vector<Real> df = Real( 0.5 ) * dt * ( Force[i] / Mass[i] );  // for testing makes velo-verlet
		 std::vector<Real> df = Real( 0.5 ) * dt * ( Force[i] / Mass[i] - NoseChi * Velo[i] );
		 Pos[i] = Pos[i] + dv + dt * df ;
		 ApplyPeriodicBound( Pos[i] );
		 std::vector<Real> v = MatVec( Box , Velo[i] );
		 Ek = Ek + DotProduct( v , v ) * Mass[i];
		 Velo[i] = Velo[i] + df ;
	 }
	 Real fs = ( Ek - NoseG * temperature ) / VirtMass;
	 NoseChi = NoseChi + Real( 0.5 ) * dt * fs ;
}




 void ThermostatsVeloVerlet::NoseHooverStepB( const std::vector<VecNSlice<Real> >& Force , 
      		                        std::vector<VecNSlice<Real> >& Velo  , 
					Real dt , 
					const std::vector<std::vector<Real> >& Box , 
											  const std::vector<Real>& Mass ){

	 Real Ek = Real( 0 );
#if USE_OPENMP
#pragma omp parallel for reduction(+:Ek)
#endif
	 for ( auto i = 0 ; i < Force.size() ; i++ ){
		 std::vector<Real> df = Real( 0.5 ) * dt * ( Force[i] / Mass[i] - NoseChi * Velo[i] );
		 //std::vector<Real> df = Real( 0.5 ) * dt * ( Force[i] / Mass[i] ); // for testing makes velo verlet
		 std::vector<Real> v = MatVec( Box , Velo[i] );
		 Ek = Ek + DotProduct( v , v ) * Mass[i];
		 Velo[i] = Velo[i] + df ;
	 }
	 Real fs = ( Ek - NoseG * temperature ) / VirtMass;
	 NoseChi = NoseChi + Real( 0.5 ) * dt * fs ;
 }
 /// end nose-hoover implemantation
 // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 


 /// Langevin implementation
 void ThermostatsVeloVerlet::LangevinStepA( std::vector<VecNSlice<Real> >& Pos , 
		                            const std::vector<VecNSlice<Real> >& Force , 
			                    std::vector<VecNSlice<Real> >& Velo , 
			                    Real dt ,
			                    const std::vector<std::vector<Real> >& Box , 
			                    const std::vector<Real>& Mass ){
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
	 for ( auto i = 0 ; i < Pos.size() ; i++ ){
		 Velo[i] = Velo[i] + Real( 0.5 ) * dt * Force[i] / Mass[i] + 
		 	             MatVec( Box , LangevinB[i] * RandomForce[i] );
		 Pos[i]  = Pos[i] + LangevinC * Velo[i];
		 ApplyPeriodicBound( Pos[i] );
     }
 }


 void ThermostatsVeloVerlet::LangevinStepB( const std::vector<VecNSlice<Real> >& Force , 
					                        std::vector<VecNSlice<Real> >& Velo , 
					                        Real dt ,
					                        const std::vector<std::vector<Real> >& Box , 
					                        const std::vector<Real>& Mass ){

#ifdef USE_OPENMP
#pragma omp parallel for
#endif
	 for ( auto i = 0 ; i < Velo.size() ; i++ ){
		 Velo[i] = LangevinA * Velo[i] + MatVec( Box , LangevinB[i] * RandomForce[i] ) + 
		 	      Real( 0.5 ) * dt * Force[i] / Mass[i];
	 }
 }


 void ThermostatsVeloVerlet::LangevinRandomNumbers( void ){
	 for ( auto i = 0 ; i < RandomForce.size() ; i++ ){
		 for ( auto j = 0 ; j < RandomForce[ i ].size() ; j++ ){
			 RandomForce[i][j]  =  Rand.generate_gauss_dist( Real( 0 ) , Real( 1 ) );
		 }
	 }
 }
//// end of langevin implementation
//// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


 // removing center of mass drift
 void ThermostatsVeloVerlet::RemoveCOMDrift( std::vector<VecNSlice<Real> >& Velo , const std::vector<Real>& Mass ){
	 std::vector<Real> COM;
	 COM.resize( Velo[0].size() );
	 Real TotMass = Real( 0 );
	 for ( auto i = 0 ; i < Velo.size() ; i++ ){
		 COM = COM + Mass[i] * Velo[i];
		 TotMass = TotMass + Mass[i];
	 }
	 COM = COM / TotMass;
	 for ( auto i = 0 ; i < Velo.size(); i++ ){
		 Velo[i] = Velo[i] - COM;
     }
 }
 //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~




 // implementation of integrators
 void ThermostatsVeloVerlet::ThermostatIntegrate( std::vector<VecNSlice<Real> >& Pos ,
                                          std::vector<VecNSlice<Real> >& Force ,
                                          std::vector<VecNSlice<Real> >& Velo ,
   	                                  Real dt ,
			                  const std::vector<std::vector<Real> >& Box , 
			                  const std::vector<std::vector<Real> >& RecBox , 
			                  ComputeForces GetForce , const std::vector<Real>& Mass ){


	 if ( ThermoType == 0 ){
		 AndersenStepA( Pos , Force , Velo , dt , Mass );   // update position
		 // replace this line to compute more forces or to compute
		 // certain ones forces.h contains special function
		 // for every application
		 //Force = GetForce.SupplyHarmonicOnly( Pos );
		 GetForce.SupplyMorseOnly( Pos , Box , RecBox , Force );
		 //Force   = GetForce.SupplyHarmonicAndOnSite( Pos );
		 AndersenStepB( Velo , Force , dt , RecBox , Mass );
	 }else if ( ThermoType == 1 ){
		 NoseHooverStepA( Pos , Force , Velo , dt , Box , Mass );
		 GetForce.SupplyMorseOnly( Pos , Box , RecBox , Force );
		 NoseHooverStepB( Force , Velo , dt , Box , Mass );
	 }else if ( ThermoType == 2 ){
		 LangevinStepA( Pos , Force , Velo , dt , RecBox , Mass );
		 GetForce.SupplyMorseOnly( Pos , Box , RecBox , Force );
		 LangevinRandomNumbers();
		 LangevinStepB( Force , Velo , dt , RecBox , Mass );
	 }else if ( ThermoType == 3 ){
		 // BDP thermostat
		 LeapFrogStepA( Pos , Velo , Force , Mass , dt );  
		 GetForce.SupplyMorseOnly( Pos , Box , RecBox , Force );
                 LeapFrogStepB( Velo , Force , Mass , dt );
		 if ( BDPHeatGrad ){
		         BDPComputeKineticEnergyHeatGrad( Velo , Box , Mass );
		         BDPRecomputeAlphaHeatGrad();
			 if ( EquiSteps <  HeatVars.EquiSteps ){
		             BDPRescaleVelocitiesHeatGrad( Velo );
			 }
			 else{
		             BDPRescaleVelocitiesHeatGradAM( Velo );
			 }
			 EquiSteps++;
		 }
		 else{
		         //thermostat 
		         BDPComputeKineticEnergy( Velo , Box , Mass );
		         BDPRecomputeAlpha();
		         BDPRescaleVelocities( Velo );
		 }
	 }

	 RemoveCOMDrift( Velo , Mass );
 }



 void ThermostatsVeloVerlet::ThermostatIntegrateHarmonicDoubleWellPerturbation( 
		 std::vector<VecNSlice<Real> >& Pos ,
                 std::vector<VecNSlice<Real> >& Force ,
                 std::vector<VecNSlice<Real> >& Velo ,
   	         Real dt ,
		 const std::vector<std::vector<Real> >& Box , 
		 const std::vector<std::vector<Real> >& RecBox , 
		 ComputeForces GetForce , const std::vector<Real>& Mass ){


	 if ( ThermoType == 0 ){
		 AndersenStepA( Pos , Force , Velo , dt , Mass );   // update position
		 // replace this line to compute more forces or to compute
		 // certain ones forces.h contains special function
		 // for every application
		 // GetForce.SupplyHarmonicOnly( Pos , Force );
		 //GetForce.SupplyDoubleWellOnly( Pos , Box , RecBox , Force );
		 //GetForce.SupplyOnSiteForceOnly( Pos , Box , RecBox , Force );
		 GetForce.SupplyHarmonicAndOnSiteForce( Pos , Box , RecBox , Force );
		 AndersenStepB( Velo , Force , dt , RecBox , Mass );
	 }else if ( ThermoType == 1 ){
		 NoseHooverStepA( Pos , Force , Velo , dt , Box , Mass );
		 // GetForce.SupplyHarmonicOnly( Pos , Force );
		 //GetForce.SupplyDoubleWellOnly( Pos , Force );
		 //GetForce.SupplyDoubleWellOnly( Pos , Box , RecBox , Force );
		 //GetForce.SupplyOnSiteForceOnly( Pos , Box , RecBox , Force );
		 GetForce.SupplyHarmonicAndOnSiteForce( Pos , Box , RecBox , Force );
		 NoseHooverStepB( Force , Velo , dt , Box , Mass );
	 }else if ( ThermoType == 2 ){
		 LangevinStepA( Pos , Force , Velo , dt , RecBox , Mass );
		 // GetForce.SupplyHarmonicOnly( Pos , Force );
		 //GetForce.SupplyDoubleWellOnly( Pos , Box , RecBox , Force );
		 //GetForce.SupplyOnSiteForceOnly( Pos , Box , RecBox , Force );
		 GetForce.SupplyHarmonicAndOnSiteForce( Pos , Box , RecBox , Force );
		 LangevinRandomNumbers();
		 LangevinStepB( Force , Velo , dt , RecBox , Mass );
	 }else if ( ThermoType == 3 ){
		 // BDP thermostat
		 LeapFrogStepA( Pos , Velo , Force , Mass , dt );  
		 //GetForce.SupplyHarmonicOnly( Pos , Box , RecBox , Force );
		 //GetForce.SupplyDoubleWellOnly( Pos , Box , RecBox , Force );
		 //GetForce.SupplyOnSiteForceOnly( Pos , Box , RecBox , Force );
		 GetForce.SupplyHarmonicAndOnSiteForce( Pos , Box , RecBox , Force );
                 LeapFrogStepB( Velo , Force , Mass , dt );
		 if ( BDPHeatGrad ){
		         BDPComputeKineticEnergyHeatGrad( Velo , Box , Mass );
		         BDPRecomputeAlphaHeatGrad();
			 if ( EquiSteps <  HeatVars.EquiSteps ){
		             BDPRescaleVelocitiesHeatGrad( Velo );
			 }
			 else{
		             BDPRescaleVelocitiesHeatGradAM( Velo );
			 }
			 EquiSteps++;
		 }
		 else{
		         //thermostat 
		         BDPComputeKineticEnergy( Velo , Box , Mass );
		         BDPRecomputeAlpha();
		         BDPRescaleVelocities( Velo );
		 }
	 }

	 RemoveCOMDrift( Velo , Mass );
 }




 void ThermostatsVeloVerlet::ThermostatIntegrateMLFF( std::vector<VecNSlice<Real> >& Pos ,
		                              std::vector<VecNSlice<Real> >& ForceDir ,
		                              std::vector<VecNSlice<Real> >& ForceCar ,
                                              std::vector<VecNSlice<Real> >& Velo ,
					      Real dt ,
					      const std::vector<std::vector<Real> >& Box , 
					      const std::vector<std::vector<Real> >& RecBox , 
					      MachineLearningForces GetForce , const std::vector<Real>& Mass ){


	 if ( ThermoType == 0 ){
		 AndersenStepA( Pos , ForceDir , Velo , dt , Mass );   // update position
		 // replace this line to compute more forces or to compute
		 // certain ones forces.h contains special function
		 // for every application
		 GetForce.ComputeDescriptors( Pos , Box );
		 GetForce.ComputeForces( ForceDir , ForceCar , RecBox );
		 AndersenStepB( Velo , ForceDir , dt , RecBox , Mass );
	 }else if ( ThermoType == 1 ){
		 NoseHooverStepA( Pos , ForceDir , Velo , dt , Box , Mass );
		 GetForce.ComputeDescriptors( Pos , Box );
		 GetForce.ComputeForces( ForceDir , ForceCar , RecBox );
		 NoseHooverStepB( ForceDir , Velo , dt , Box , Mass );
	 }else if ( ThermoType == 2 ){
		 LangevinStepA( Pos , ForceDir , Velo , dt , RecBox , Mass );
		 GetForce.ComputeDescriptors( Pos , Box );
		 GetForce.ComputeForces( ForceDir , ForceCar , RecBox );
		 LangevinRandomNumbers();
		 LangevinStepB( ForceDir , Velo , dt , RecBox , Mass );
	 }

	 RemoveCOMDrift( Velo , Mass );
 }

 // end integrators
 // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~










 void ThermostatsVeloVerlet::AddParticle( Real value ){
	 /*if ( ThermoType == 0 ){

	 }else if ( ThermoType == 1 ){

	 }else if ( ThermoType == 2 ){
		 LangevinB.push_back( value );
	     std::vector<Real> data;
		 for ( auto i = 0 ; i < RandomForceData[0].size() ; i++ ){
			 data.push_back( Rand.generate_gauss_dist( Real( 0 ) , Real( 1 ) ) );
		 }
		 RandomForceData.push_back( data );
	 }*/
 }
