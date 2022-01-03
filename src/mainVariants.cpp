#include "mainVariants.h"

void MainStyles::SetUpSimulation( AtomType & atoms , SystemParams & System , std::string inputFile){

	Units unit;

	Ran Rand( std::chrono::system_clock::to_time_t( std::chrono::system_clock::now() ) );
	
	FlagFinder Find( inputFile );

	System.dimension  =  int( Find.CheckFlag( "Ndim" ) );

	std::cout <<"##################################"<< std::endl;
	std::cout <<"########## PARAMETERS ############"<< std::endl;
	std::cout << "dimension " << System.dimension << std::endl;

	std::vector<Real> temp =  Find.CheckFlagVec( "NCells" );
	System.NParticles =  std::vector<int> ( temp.begin() , temp.end() );

	std::cout << "Repeating units " << System.NParticles[0] << "   "
		                        << System.NParticles[1] << "   "
                                        << System.NParticles[2] << std::endl;

	System.deltaX = Find.CheckFlagVec( "DeltaX" );
	std::cout << "Grid spacing [\\AA]" << std::setw( 8 ) << System.deltaX[0] << "   "
	   	                           << std::setw( 8 ) << System.deltaX[1] << "   "
                                           << std::setw( 8 ) << System.deltaX[2] << std::endl;


	System.Box = std::vector<std::vector<Real> > ( System.dimension , std::vector<Real> ( System.dimension , Real( 0.0 ) ) );
	System.RecBox = std::vector<std::vector<Real> > ( System.dimension , std::vector<Real> ( System.dimension , Real( 0.0 ) ) );

	System.timeStep = Find.CheckFlag( "TStep" ) * unit.timeConv;
	std::cout << "Time step in ps " << System.timeStep / unit.timeConv << std::endl;

	System.temperature = Find.CheckFlag( "Temp" );
	std::cout << "Temperature [K] " << System.temperature << std::endl;
	

        // determine cutoff radius in an Angstroem	
	System.Cutoff = Find.CheckFlag( "Cutoff" );
	std::cout << "Cutoff distance set to " << System.Cutoff << "[\AA]" << std::endl;


	std::cout << "Analysis frequency is set to" << std::endl;
	System.AnaFrequ =  int ( Find.CheckFlag( "AnaFrequ" ) );

	for ( auto i = 0 ; i < System.deltaX.size() ; i++ ){
		System.Box[ i ][ i ] = Real( System.NParticles[i] ) *
			               System.deltaX[i];
		System.RecBox[ i ][ i ] = Real( 1 ) / System.Box[ i ][ i ];
        }


	System.Natoms = std::accumulate( std::begin( System.NParticles ) ,
			                 std::end( System.NParticles ) , 1 ,
					 std::multiplies<Real>() );

	System.Nsteps = int( Find.CheckFlag( "NSTEPS" ) );

	std::cout << "Number of MD steps " << System.Nsteps << std::endl;


	atoms.PosDirData.resize( System.Natoms * System.dimension );
	atoms.InitializePointers( atoms.PosDirData , atoms.PosDir , System.dimension );

	atoms.PosCarData.resize( System.Natoms * System.dimension );
	atoms.InitializePointers( atoms.PosCarData , atoms.PosCar , System.dimension );

	atoms.VeloDirData.resize( System.Natoms * System.dimension );
	atoms.InitializePointers( atoms.VeloDirData , atoms.VeloDir , System.dimension );

	atoms.VeloCarData.resize( System.Natoms * System.dimension );
	atoms.InitializePointers( atoms.VeloCarData , atoms.VeloCar , System.dimension );

	atoms.ForceDirData.resize( System.Natoms * System.dimension );
	atoms.InitializePointers( atoms.ForceDirData , atoms.ForceDir , System.dimension );
	
	atoms.ForceCarData.resize( System.Natoms * System.dimension );
	atoms.InitializePointers( atoms.ForceCarData , atoms.ForceCar , System.dimension );

	atoms.Mass.resize( System.Natoms );
	atoms.Type.resize( System.Natoms );
	atoms.Spec.resize( System.Natoms );
	atoms.Epot.resize( System.Natoms );
	atoms.Ekin.resize( System.Natoms );

	// reading thermostat stuff
	System.Thermostyle  =  int( Find.CheckFlag( "ThermST" ) );
	if ( System.Thermostyle == 0 ){
		std::cout << "Andersen thermostat is used" << std::endl;
		System.AndersenFrequ  =  Find.CheckFlag( "COLFA" ) / unit.timeConv ;
		std::cout << "Andersen collision frequ " << std::setw( 8 ) <<
		System.AndersenFrequ * System.timeStep <<std::endl;
        }
	else if ( System.Thermostyle == 1 ){
		std::cout << "Nose-Hoover thermostat is used" << std::endl;
		System.NoseMass  =  Find.CheckFlag( "NoseM" ) * unit.timeConv ;
		std::cout << "Nose-Hoover mass " << std::setw( 8 ) <<
		              System.NoseMass << std::endl;
        } 
	else if ( System.Thermostyle == 2 ){
		std::cout << "Langevin thermostat is used" << std::endl;
		System.LangevinFriction  =  Find.CheckFlag( "LangFric" ) * unit.timeConv ;
		std::cout << "Langevin friction coefficient " << std::setw( 8 ) <<
		              System.LangevinFriction << std::endl;
        }
	else if ( System.Thermostyle == 3 ){
		std::cout << "Bussi-Donadio-Parrinello thermostat is used" << std::endl;
		System.BDPThermoTau  =  Find.CheckFlag( "BDPTau" ) * unit.timeConv ;
		std::cout << "BDP time scale tau " << std::setw( 8 ) <<
		              System.BDPThermoTau / unit.timeConv << std::endl;
	}




        // delta vector for phonon dispersion in
	// harmonic approximation
	// supply in cartesian coordinates	
	System.DeltaPhonon = Find.CheckFlagVec( "DeltaPhon" );
	if ( System.DeltaPhonon.size() == 0 ){
		System.PhononOnOff = false;
        }
	else{
		System.PhononOnOff = true;
		System.DeltaPhonon = MatVec( System.RecBox , System.DeltaPhonon );
	}


	System.Output  =  0;
	System.Output  =  int( Find.CheckFlag( "Output" ) );
	if ( System.Output < 0 ){
		System.Output = 0;
        }


	Real M = Find.CheckFlag( "mass" );
	std::cout << "Particle masses [Da] " << std::setw( 8 ) <<
		         M <<std::endl;


	System.StartSample = int( Find.CheckFlag( "EquiStep" )  );

	if ( System.StartSample > 0 ){
		std::cout <<  "Starting to take samples at MD step "
			      <<  System.StartSample << std::endl;
        }else{
		std::cout << "No equibrilation time set; start sampling at first step" << std::endl;
		System.StartSample = -1;
        }

	System.SwitchThermo = int( Find.CheckFlag( "SwitchNVE" ) );
	if ( System.SwitchThermo > 0 ){
	      std::cout << "Switching to NVE ensemble after " << 
	      System.SwitchThermo << " steps activated" << std::endl;
	}


	System.NNMorseParams = Find.CheckFlagVec( "NNMorse" );
	if ( System.NNMorseParams.size() != 3 ){
		System.NNMorseParams.resize( 3 );
	        System.NNMorseParams[0] = 0.2348 ; // depth
	        System.NNMorseParams[1] = 1.1836;  // steepness
	        System.NNMorseParams[2] = 3.733 ;  // equilibrium distance
	}
	std::cout << "Nearest neighbour Morse potential parameters" << std::endl;
        std::cout <<  System.NNMorseParams[ 0 ]  <<  "    "
	          <<  System.NNMorseParams[ 1 ]  <<  "    " 
	          <<  System.NNMorseParams[ 2 ]  <<  "    " 
	          <<  std::endl;
	std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;


	System.RestartFile =  Find.CheckStringFlag( "restart_file" );
	unsigned int len   =  System.RestartFile.length();

	if ( !System.RestartFile.substr( 0, len ).compare( "NOTFOUND" ) ){
	    unsigned int col = 0;
	    for ( auto i = 0 ; i < System.NParticles[0] ; i++ ){
	      	Real x = Real( i ) * System.deltaX[0] + System.deltaX[0] * Real( 0.5 );
	       	for ( auto j = 0 ; j < System.NParticles[1] ; j++ ){
	       	    Real y = Real( j ) * System.deltaX[1] + System.deltaX[1] * Real( 0.5 );
	       		for ( auto k = 0 ; k < System.NParticles[2] ; k++ ){
	       			Real z = Real( k ) * System.deltaX[2] + System.deltaX[2] * Real( 0.5 );
	       			atoms.Type[ col ] = 1;
		   		    if ( !System.PhononOnOff ){
					    //atoms.PosCar[ col ][ 0 ] = x + 0.1 * ( Rand.doub() - 1.0 )*System.deltaX[0];
	       			            //atoms.PosCar[ col ][ 1 ] = y + 0.1 * ( Rand.doub() - 1.0 )*System.deltaX[1];
	       			            //atoms.PosCar[ col ][ 2 ] = z + 0.1 * ( Rand.doub() - 1.0 )*System.deltaX[2];
	       			            atoms.PosCar[ col ][ 0 ] = x ;
	       			            atoms.PosCar[ col ][ 1 ] = y ;
	       			            atoms.PosCar[ col ][ 2 ] = z ;
	       			            //atoms.PosCar[ col ][ 0 ] = x + 0.5 * Rand.generate_gauss_dist( float( 0 ) , 
				            //     		                                             0.5 * System.temperature * unit.kB );
	       			            //atoms.PosCar[ col ][ 1 ] = y + 0.5 * Rand.generate_gauss_dist( float( 0 ) , 
				            //     		                                             0.5 * System.temperature * unit.kB );
	       			            //atoms.PosCar[ col ][ 2 ] = z + 0.5 * Rand.generate_gauss_dist( float( 0 ) ,
				            //		                                              0.5 * System.temperature * unit.kB );	   
		   		    }
		   		    else{
					    atoms.PosCar[ col ][ 0 ] = x ;
	       			            atoms.PosCar[ col ][ 1 ] = y ;
	       			            atoms.PosCar[ col ][ 2 ] = z ;
		   		   }
	       			   atoms.PosDir[ col ] = MatVec( System.RecBox , atoms.PosCar[ col ] );
	       			   atoms.Mass[col] = M;
	       			   col++;
	       	    }
	       	}
	    }
        }else{
		  atoms.ReadRestartStructure( System.RestartFile );
		  std::cout << "Old structure read " << std::endl;
		  std::cout << "Check that your other parameters are in" << std::endl;
		  std::cout << "agreement with the supplied structure" << std::endl;
		  for ( auto i = 0 ; i < atoms.PosDir.size() ; i++ ){
			  std::cout << atoms.PosDir[i][0] << "   "
		  	            << atoms.PosDir[i][1] << "    "
	   	     	    	    << atoms.PosDir[i][2] << std::endl;
	          atoms.Mass[ i ] = M;
		  }
	}


	System.ParticleInsert = ( int ) ( Find.CheckFlag( "ParticleInsert" )  );
	if ( System.ParticleInsert > 0 ){
		std::cout << "New particle will be inserted every "
			       << System.ParticleInsert << " steps" << std::endl;

	}else{
		std::cout << "Particle insertion switched off" << std::endl;
	}

	
	// read heat transport variables 
	System.HeatVars.HeatGradient = ( int ) ( Find.CheckFlag( "HeatGrad" )  );
	if ( System.HeatVars.HeatGradient > 0 ){
		std::cout << "Heat gradient was switched on" << std::endl;
	        System.HeatVars.temperatureA =  Find.CheckFlag( "TempA" );
	        System.HeatVars.temperatureB =  Find.CheckFlag( "TempB" );
		if ( System.HeatVars.temperatureA < 0 ){
			throw ( System.HeatVars.temperatureA );
		}
		std::cout << "Temperature A for heat transport set to" 
			  <<  std::setw( 8 ) << System.HeatVars.temperatureA << std::endl;
		if ( System.HeatVars.temperatureB < 0 ){
			throw ( System.HeatVars.temperatureB );
		}
		std::cout << "Temperature B for heat transport set to" 
			  <<  std::setw( 8 ) << System.HeatVars.temperatureB << std::endl;

	        System.HeatVars.Direction    =  ( int ) ( Find.CheckFlag( "HeatDir" ) );
		if ( System.HeatVars.Direction < 1 or System.HeatVars.Direction > 3 ){
			std::cout << "Direction for temperature gradient not given" << std::endl;
			std::cout << "Supply a number {1,2,3}" << std::endl;
			std::cout << "x ---> 1" << std::endl;
			std::cout << "y ---> 2" << std::endl;
			std::cout << "z ---> 3" << std::endl;
		}
		else{
			std::cout << "HeatGardient direction set to " 
				  << System.HeatVars.Direction << std::endl;
		}

		System.HeatVars.ActiveMeasure = ( int ) ( Find.CheckFlag( "ActMeasure" ) );
		if ( System.HeatVars.ActiveMeasure < 0 ){
			std::cout << "Active heat gradient measure switched off" << std::endl;
		}
		else{
			std::cout << "Active heat gradient measure switched on" << std::endl;
		}

		System.HeatVars.EquiSteps  =  ( int ) ( Find.CheckFlag( "ActMEqui" ) );
	}

	std::cout <<"##################################"<< std::endl;
}




int MainStyles::MorseSolidMain( const char * fname , const unsigned int SimType ){

	AtomType atoms;       // contains information about atoms
	SystemParams System;  // information about the simulation conditions

	SetUpSimulation( atoms , System , fname );


	if ( SimType == 2 ){
		// allows for liquid and melting simulations
		return MorseListUpdate( atoms , System );
	}






	// set up force class
	ComputeForces Force( atoms.PosDir , ( unsigned int )( 6 ) , System.Box );
	// set up harmonic force  coupling strength equilibrium dist cart equilibrium dist dir
	Force.SetUpHarmonicPart( Real( 5.000 ) , System.Box[0][0] / Real( System.NParticles[0] ) );
	// set up morse potential 
        Force.SetUpMorsePart( System.NNMorseParams[0] , System.NNMorseParams[1] , System.NNMorseParams[2] );


        ThermostatsVeloVerlet Integrator;
	// set up thermostat integrator velocity verlet
	if ( System.HeatVars.HeatGradient < 1 ){
	    if ( System.Thermostyle == 0 ){
	            Integrator.InitAndersen( System.temperature ,
				System.AndersenFrequ*System.timeStep , atoms.Mass );
	    }
	    else if ( System.Thermostyle == 1 ){
	    	Integrator.InitNoseHoover( System.NoseMass , 
	          	   atoms.PosDir.size() * atoms.PosDir[0].size() + 1 ,
			   System.temperature , atoms.Mass );
	    }
	    else if ( System.Thermostyle == 2 ){
	    	Integrator.InitLangevin( System.LangevinFriction ,
				         System.temperature , 
	    			         System.timeStep , System.dimension , atoms.Mass );
	    }
	    else if ( System.Thermostyle == 3 ){
		    Integrator.InitBDP( System.BDPThermoTau , System.temperature ,
				        System.timeStep , 
					atoms.PosDir.size() * atoms.PosDir[0].size() , // total degrees of freedom
					atoms.Mass );
	    }
	}
	// do the heat gradient initialization
	else{
	    if ( System.Thermostyle == 0 ){
		System.HeatVars.nu = System.AndersenFrequ * System.timeStep;
		Integrator.InitHeatGradientThermoAndersen( System.HeatVars , 
				System.NParticles , atoms.Mass );
	    }
	    else if ( System.Thermostyle == 2 ){
	    	Integrator.InitHeatGradientThermoLangevin( 
				         System.HeatVars , System.LangevinFriction ,
	    			         System.timeStep , System.dimension , 
					 System.NParticles , atoms.Mass );
	    }
	    else if ( System.Thermostyle == 3 ){
		    Integrator.InitHeatGradientThermoBDP( System.HeatVars , System.BDPThermoTau ,
				System.timeStep , atoms.Mass , System.NParticles , 
				System.NParticles.size() );
	    }
	}

	// force initialization function call has to match
	// the one used for force computation in andersen routine
	//atoms.ForceDir =  Force.SupplyHarmonicOnly( atoms.PosDir );
	Force.SupplyMorseOnly( atoms.PosDir , System.Box , System.RecBox , atoms.ForceDir );
	//atoms.ForceDir = Force.SupplyHarmonicAndOnSite( atoms.PosDir );



	//phonon finite differences is computed here
	if ( System.PhononOnOff ){
		FDPhonons Phonons( System.Box , System.RecBox , System.NParticles , System.DeltaPhonon );
		Phonons.main( atoms.PosDir , Force , System.Box , atoms.Mass );
		return 0;
	}



	std::vector<std::string> types;
	types = { "H" };
	std::vector<int> Nats ( 1 , System.Natoms );

	// STRUCTURE output
	FILE * xdatcar;
	if ( System.Output == 0 ){
		xdatcar = fopen( "XDATCAR" , "w" );
	}else if ( System.Output == 1 ){
		xdatcar = fopen( "traj.lammpstrj" , "w" );
        }
	else if ( System.Output == 2 ){
		xdatcar = fopen( "md.out" , "w" );
	}

	//FILE * ForceCar;
	//ForceCar = fopen( "Forces" , "w" );


        WriteOutput StrucOut( xdatcar );
        //WriteOutput ForceOut( ForceCar );
	// xdat init
	//StrucOut.XdatInit( System.Box , types , Nats );
	//StrucOut.AddStruc( atoms.PosDir );

	if ( System.Output == 0 ){
		StrucOut.XdatInit( System.Box , types , Nats );
	        StrucOut.AddStruc( atoms.PosDir );
        }else if ( System.Output == 1 ){
            	StrucOut.LammpsWrite( atoms.PosDir , System.Box , atoms.Type ,
            	               System.Natoms );
        }
            else if ( System.Output == 2 ){
            	StrucOut.PhQWriteStruc( atoms.PosDir , true , System.Nsteps );
        }

	// parameter output
	FILE * Efile;  // energy file
	Efile = fopen( "energy.out" , "w" );
	WriteOutput MDSpec( Efile );

	std::vector<Real> MDData;
	MDData.resize( 4 );


	//velocity init
	Integrator.InitVelocities( atoms.VeloCar , System.dimension , System.Box , atoms.Mass );
	atoms.Transform( atoms.VeloDir, atoms.VeloCar , System.RecBox );



        //System.Epot  =  Force.PotentialEnergyHarmonic( atoms.PosDir , System.Box );
        System.Epot  =  Force.PotentialEnergyMorse( atoms.PosDir , System.Box );
        //System.Epot  =  Force.PotentialEnergyHarmonicAndOnSite( atoms.PosDir , System.Box );
	System.Ekin  =  ComputeKineticEnergy( atoms.VeloCar , atoms.Mass , atoms.Ekin );
	System.Etot  =  System.Epot + System.Ekin ;
	System.Tact  =  ComputeTemperature( System.Ekin , System.Natoms );

	MDData = { System.Tact , System.Ekin , System.Epot , System.Etot };
	//MDSpec.WriteVector( MDData );


	// initialize analysis
	PairDistFunction DistDistributions;
	DistDistributions.DistCenter( System.Box[0][0] / Real( 6 ) , atoms.PosDir , 250 );
	DistDistributions.PairDistInit( System.Box[0][0] , 250 );
	AngularDist3d AngDist( 50 , 100 );
	

	AnalyzePotential PotentialAna( atoms.PosDir );


        std::vector<Real> COMVel;
        COMVel.resize( atoms.VeloDir[0].size() );
	// this part can be used to test the potential energy computation

//	{
//
//	   int Steps    =  100;
//	   Real dstep   =  atoms.PosDir[1][2] - atoms.PosDir[0][2];
//	   Real EquiPos =  atoms.PosDir[0][2] * System.Box[0][0];
//           dstep        =  dstep  /  Real( Steps );
//
//	   atoms.PosDir[0][2]   =  atoms.PosDir[0][2] - dstep * Real( Steps ) / Real( 2 );
//
//
//	   FILE * Efile;  // energy file
//	   Efile = fopen( "test.out" , "w" );
//
//
//	   for ( auto i = 0 ; i < Steps ; i++ ){
//	   	System.Epot    =  Force.PotentialEnergyMorse( atoms.PosDir , System.Box );
//	   	Force.SupplyMorseOnly( atoms.PosDir , System.Box , System.RecBox , atoms.ForceDir );
//		fprintf( Efile , "   %f" , atoms.PosDir[0][2] * System.Box[0][0] - EquiPos );
//		fprintf( Efile , "   %f" , System.Epot );
//		fprintf( Efile , "   %f" , atoms.ForceDir[0][2] * System.Box[0][0] );
//		fprintf( Efile , "\n" );
//                fflush( Efile );
//	   	atoms.PosDir[0][2]  =  atoms.PosDir[0][2] + dstep;
//	   }
//
//	throw;
//	}

	for ( auto i = 0 ; i < System.Nsteps ; i++ ){


		//Integrator.AndersenIntegrate( atoms.PosDir , atoms.ForceDir , atoms.VeloDir ,
		//            System.timeStep , MexPot , HarmonicCoupling , LinkLister.NNList , System.RecBox );
		//        Calling anderson integrator to do single step
		//Integrator.VelocityExchange( atoms.VeloDir );
		Integrator.ThermostatIntegrate( atoms.PosDir , atoms.ForceDir , atoms.VeloDir ,
	                         	        System.timeStep , System.Box , 
						System.RecBox , Force , atoms.Mass );

		if ( System.SwitchThermo == i ){
			std::cout << "Switched to NVE ensemble at step " << i << std::endl;
			System.Thermostyle = 0;
			System.AndersenFrequ = Real( 0.0 );
          	        Integrator.InitAndersen( System.temperature , 
					System.AndersenFrequ*System.timeStep , atoms.Mass );
	        }

		// analysis and output
		if ( i%System.AnaFrequ == 0  and i > System.StartSample ){
		   //if nvt is used to equibrilate one can switch to nve after N steps

		   printf( "Evaluating sample Nr , %i \n" , i );
	           // compute pair dist
		   DistDistributions.ComputeDistToCenter( atoms.PosDir , System.Box );
		   //DistDistributions.ComputePairDistribution( atoms.PosDir , System.Box );
		   AngDist.ComputeAngularCOMCentered( atoms.PosDir , DistDistributions.Center , System.Box );

	           atoms.Transform( atoms.VeloCar, atoms.VeloDir , System.Box );
	           //atoms.Transform( atoms.PosCar, atoms.PosDir , System.Box );
		   //std::vector<Real> AvVel;
		   //AvVel.resize( atoms.VeloDir[0].size() );
		   //for ( auto a = 0 ; a < atoms.VeloDir.size() ; a++ ){
		   //	   AvVel = AvVel + atoms.VeloDir[a];
		   //}
		   //AvVel = MatVec( System.Box , AvVel ) / Real( atoms.VeloDir.size() );
		   //COMVel = COMVel + AvVel * System.timeStep;

	            //atoms.Transform( atoms.ForceCar, atoms.ForceDir , System.Box );
		    //ForceOut.AddStruc( atoms.ForceCar );
		    //
		    //
		    //potential analysis

		    Force.ComputeMorsePotentialEnergyArray( atoms.PosDir , System.Box , atoms.Epot );
		    //PotentialAna.PotentialAnaMain( atoms.PosDir , atoms.Epot , System.Box , atoms.ForceDir );



		    // compute md specs
		    //System.Epot   =  Force.PotentialEnergyHarmonic( atoms.PosDir , System.Box );
		    System.Epot  =  Force.PotentialEnergyMorse( atoms.PosDir , System.Box );
		    //System.Epot  =  Force.PotentialEnergyHarmonicAndOnSite( atoms.PosDir , System.Box );
	            System.Ekin  =  ComputeKineticEnergy( atoms.VeloCar , atoms.Mass , atoms.Ekin );
	            System.Etot  =  System.Epot + System.Ekin ;
	            System.Tact  =  ComputeTemperature( System.Ekin , System.Natoms );
		    
		    // compute layer temperatures and write them to thermostats
		    // debug file
		    
		    //Integrator.ComputeLayerTemperature( atoms.VeloCar , atoms.Mass );

		    // write md specs
	            //MDData = {System.Tact , System.Ekin , System.Epot, System.Etot , COMVel[0],COMVel[1],COMVel[2] };
	            MDData = {System.Tact , System.Ekin , System.Epot, System.Etot };
	            MDSpec.WriteVector( MDData );
		    // write structure
		    if ( System.Output == 0 ){
		            StrucOut.AddStruc( atoms.PosDir );    //XDAT format
		    }
		    //Lammps
		    else if ( System.Output == 1 ){
		            StrucOut.LammpsWrite( atoms.PosDir , System.Box , atoms.Type ,
		                          System.Natoms );
		    }
		    if ( System.Output == 2 ){
		            StrucOut.PhQWriteStruc( atoms.PosDir , false , i );
		    }
		}
        }
	//DistDistributions.WriteOutputCOM( "CenterDist.out" );
	//DistDistributions.WriteOutputPD( "PairDist.out" );
	//AngDist.WriteOutput( "PolarDist.out" );
	//PotentialAna.PotentialAnaFinalize();

	
    return 0;	
}


  /******************************************************
     Morse potential with list update allows for melting
     of box
  ******************************************************/

 int MainStyles::MorseListUpdate( AtomType & atoms , SystemParams & System ){
	
	 
	std::vector<std::string> types;
	types = { "H" };
	std::vector<int> Nats ( 1 , System.Natoms );
	Units units;

	
	std::cout << "Morse potential with list update" << std::endl;
	ComputeForces Force( atoms.PosCar , atoms.PosDir , System.Cutoff , System.Box );


	
        Force.SetUpMorsePart( System.NNMorseParams[0] , System.NNMorseParams[1] , System.NNMorseParams[2] );
    
	// integrator
	ThermostatsVeloVerlet Integrator;
	// set up thermostat integrator velocity verlet
	if ( System.Thermostyle == 0 ){
  	        Integrator.InitAndersen( System.temperature , System.AndersenFrequ*System.timeStep , atoms.Mass );
	}
	else if ( System.Thermostyle == 1 ){
		Integrator.InitNoseHoover( System.NoseMass , 
	      	   atoms.PosDir.size() * atoms.PosDir[0].size() + 1 , System.temperature , atoms.Mass );
	}
	else if ( System.Thermostyle == 2 ){
		Integrator.InitLangevin( System.LangevinFriction , System.temperature , 
				                 System.timeStep , System.dimension , atoms.Mass );
	}
	else if ( System.Thermostyle == 3 ){
		Integrator.InitBDP( System.BDPThermoTau , System.temperature ,
		                    System.timeStep , 
				    atoms.PosDir.size() * atoms.PosDir[0].size() , // total degrees of freedom
				    atoms.Mass );
	}


	// force initialization function call has to match
	// the one used for force computation in andersen routine
	//atoms.ForceDir =  Force.SupplyHarmonicOnly( atoms.PosDir );
	Force.SupplyMorseOnly( atoms.PosDir , System.Box , System.RecBox , atoms.ForceDir );


	// STRUCTURE output
	FILE * xdatcar;
	if ( System.Output == 0 ){
		xdatcar = fopen( "XDATCAR" , "w" );
	}else if ( System.Output == 1 ){
		xdatcar = fopen( "traj.lammpstrj" , "w" );
        }
	else if ( System.Output == 2 ){
		xdatcar = fopen( "md.out" , "w" );
	}

	//FILE * ForceCar;
	//ForceCar = fopen( "Forces" , "w" );


        WriteOutput StrucOut( xdatcar );
        //WriteOutput ForceOut( ForceCar );
	// xdat init
	//StrucOut.XdatInit( System.Box , types , Nats );
	//StrucOut.AddStruc( atoms.PosDir );

	if ( System.Output == 0 ){
		StrucOut.XdatInit( System.Box , types , Nats );
	        StrucOut.AddStruc( atoms.PosDir );
        }else if ( System.Output == 1 ){
            	StrucOut.LammpsWrite( atoms.PosDir , System.Box , atoms.Type ,
            	               System.Natoms );
        }
            else if ( System.Output == 2 ){
            	StrucOut.PhQWriteStruc( atoms.PosDir , true , System.Nsteps );
        }

	// parameter output
	FILE * Efile;  // energy file
	Efile = fopen( "energy.out" , "w" );
	WriteOutput MDSpec( Efile );

	std::vector<Real> MDData;
	MDData.resize( 4 );
	
	
	//velocity init
	Integrator.InitVelocities( atoms.VeloCar , System.dimension , System.Box , atoms.Mass );
	atoms.Transform( atoms.VeloDir, atoms.VeloCar , System.RecBox );

        //System.Epot  =  Force.PotentialEnergyHarmonic( atoms.PosDir , System.Box );
        System.Epot  =  Force.PotentialEnergyMorse( atoms.PosDir , System.Box );
        //System.Epot  =  Force.PotentialEnergyHarmonicAndOnSite( atoms.PosDir , System.Box );
	System.Ekin  =  ComputeKineticEnergy( atoms.VeloCar , atoms.Mass , atoms.Ekin );
	System.Etot  =  System.Epot + System.Ekin ;
	System.Tact  =  ComputeTemperature( System.Ekin , System.Natoms );

	MDData = { System.Tact , System.Ekin , System.Epot , System.Etot };
	//MDSpec.WriteVector( MDData );


	// initialize analysis
	PairDistFunction DistDistributions;
	DistDistributions.DistCenter( System.Box[0][0] / Real( 6 ) , atoms.PosDir , 250 );
	DistDistributions.PairDistInit( System.Box[0][0] , 250 );
	AngularDist3d AngDist( 50 , 100 );


        std::vector<Real> COMVel;
        COMVel.resize( atoms.VeloDir[0].size() );


	InsertParticle Insertion( 100 , System.temperature );


	float t1 = omp_get_wtime();
	for ( auto i = 0 ; i < System.Nsteps ; i++ ){
	
		Integrator.ThermostatIntegrate( atoms.PosDir , atoms.ForceDir , atoms.VeloDir ,
	                         	        System.timeStep , System.Box , System.RecBox , Force , atoms.Mass );

		// switch to NVE ensemble
		if ( System.SwitchThermo == i ){
			std::cout << "Switched to NVE ensemble at step " << i << std::endl;
			System.Thermostyle = 0;
			System.AndersenFrequ = Real( 0.0 );
	                Integrator.InitAndersen( System.temperature , 
					System.AndersenFrequ*System.timeStep , atoms.Mass );
	         }


		// particle insertion
		//if ( i%System.ParticleInsert == 0 ){
		//	std::cout << "Insert particle at step " << i << std::endl;
		//	std::vector<Real> NewForce = Insertion.MonteCarloInsertion( atoms.PosDir , atoms.VeloCar , System.Cutoff , System.Box , Force , 
		//			                       sqrt( System.temperature * units.kB / atoms.Mass[0] ) );
		//	atoms.VeloDir.push_back( MatVec( System.RecBox , atoms.VeloCar.back() ) );
		//	atoms.Mass.push_back( atoms.Mass[ 0 ] );
		//	atoms.ForceDir.push_back( NewForce );
		//	atoms.Type.push_back( 1 );
		//    Integrator.AddParticle( System.LangevinFriction );	
		//	std::cout << "particle inserted" << std::endl;
		//}

		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ neighbor list update ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	        atoms.Transform( atoms.VeloCar, atoms.VeloDir , System.Box );
	        atoms.Transform( atoms.PosCar, atoms.PosDir , System.Box );
		Force.LinkList.UpdateNNVerletList( atoms.PosDir , atoms.PosCar , atoms.VeloCar ,
			                               System.timeStep , System.Cutoff , System.Box );
        //Force.LinkList.ComputeNNListPeriodicBound( atoms.PosDir , System.Box , System.Cutoff );


		// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

		// analysis and output
		if ( i%System.AnaFrequ == 0  and i > System.StartSample ){
			//if nvt is used to equibrilate one can switch to nve after N steps

		        printf( "Evaluating sample Nr , %i \n" , i );
	                // compute pair dist
		        DistDistributions.ComputeDistToCenter( atoms.PosDir , System.Box );
		        AngDist.ComputeAngularCOMCentered( atoms.PosDir ,
					DistDistributions.Center , System.Box );


	                atoms.Transform( atoms.ForceCar, atoms.ForceDir , System.Box );
		        //// compute md specs
		        System.Epot     =  Force.PotentialEnergyMorse( atoms.PosDir , System.Box );
	                System.Ekin     =  ComputeKineticEnergy( atoms.VeloCar , atoms.Mass , atoms.Ekin );
	                System.Etot     =  System.Epot + System.Ekin ;
	                System.Tact     =  ComputeTemperature( System.Ekin , System.Natoms );
		        System.Pressure =  ComputePressure( atoms.PosDir , atoms.ForceCar ,
		         	                               System.Box , System.temperature,
		             			       Force.LinkList.NNList );

		        //// write md specs
	                MDData = { System.Tact , System.Ekin , System.Pressure , System.Epot, System.Etot };
	                MDSpec.WriteVector( MDData );
		        // write structure
		        if ( System.Output == 0 ){
		            StrucOut.AddStruc( atoms.PosDir );    //XDAT format
		        }
		        //Lammps
		        else if ( System.Output == 1 ){
		            StrucOut.LammpsWrite( atoms.PosDir , System.Box , atoms.Type ,
		                                  System.Natoms );
		        }
		        if ( System.Output == 2 ){
		            StrucOut.PhQWriteStruc( atoms.PosDir , false , i );
		        }
		}
    }

	float t2 = omp_get_wtime();

	std::cout << "Execution time " << t2 - t1 << std::endl;

	DistDistributions.WriteOutputCOM( "CenterDist.out" );
	//DistDistributions.WriteOutputPD( "PairDist.out" );
	AngDist.WriteOutput( "PolarDist.out" );

    return 0;	

 }





 int MainStyles::MachineLearningFF( const char * fname ){
	
	 AtomType atoms;       // contains information about atoms
	 SystemParams System;  // information about the simulation conditions 
	 SetUpSimulation( atoms , System , fname );




	 std::vector<std::string> types;
 	 types = { "H" };
	 std::vector<int> Nats ( 1 , System.Natoms );
	 Units units;

	 FlagFinder Find( fname );
	 std::string WeightsFile = Find.CheckStringFlag( "WF" );
	 std::string DescriptorFile = Find.CheckStringFlag( "DF" );
	 std::string DescriptorDerivativeFile = Find.CheckStringFlag( "DDF" );

	 std::string StartFile = Find.CheckStringFlag( "Start" );

	 if ( StartFile != "NOTFOUND" ){
		 atoms.ReadStartStruc( StartFile , System );
	 }


	 // read parameters for machine learning force field
	 using namespace std::placeholders;
	 auto fp = std::bind( &AtomType::InitializePointers , atoms , _1 , _2 , _3 );
	 MachineLearningForces MLFF( WeightsFile , DescriptorFile ,
			                      DescriptorDerivativeFile , fp , atoms.PosDir.size() , System.Cutoff, 
								  atoms.PosDir , atoms.PosCar , System.Box );




	 MLFF.ComputeDescriptors( atoms.PosDir , System.Box );
	 MLFF.ComputeForces( atoms.ForceDir , atoms.ForceCar , System.RecBox );


	 // integrator
	 ThermostatsVeloVerlet Integrator;
	 // set up thermostat integrator velocity verlet
	 if ( System.Thermostyle == 0 ){
  	         Integrator.InitAndersen( System.temperature , System.AndersenFrequ*System.timeStep , atoms.Mass );
	 }
	 else if ( System.Thermostyle == 1 ){
	 	//ThermostatsVeloVerlet Integrator( System.NoseMass , atoms.PosDir.size() * atoms.PosDir[0].size() + 1 , System.temperature );
	 	Integrator.InitNoseHoover( System.NoseMass , 
	       	   atoms.PosDir.size() * atoms.PosDir[0].size() + 1 , System.temperature , atoms.Mass );
	 }
	 else if ( System.Thermostyle == 2 ){
	 	Integrator.InitLangevin( System.LangevinFriction , System.temperature , 
	 			                 System.timeStep , System.dimension , atoms.Mass );
	 }


	 // STRUCTURE output
	 FILE * xdatcar;
	 if ( System.Output == 0 ){
	 	xdatcar = fopen( "XDATCAR" , "w" );
	 }else if ( System.Output == 1 ){
	 	xdatcar = fopen( "traj.lammpstrj" , "w" );
     }
	 else if ( System.Output == 2 ){
	 	xdatcar = fopen( "md.out" , "w" );
	 }


     WriteOutput StrucOut( xdatcar );
     //WriteOutput ForceOut( ForceCar );
	 // xdat init
	 //StrucOut.XdatInit( System.Box , types , Nats );
	 //StrucOut.AddStruc( atoms.PosDir );

	 if ( System.Output == 0 ){
	 	StrucOut.XdatInit( System.Box , types , Nats );
	         StrucOut.AddStruc( atoms.PosDir );
         }else if ( System.Output == 1 ){
	 	StrucOut.LammpsWrite( atoms.PosDir , System.Box , atoms.Type ,
	 	               System.Natoms );
         }
	 else if ( System.Output == 2 ){
	 	StrucOut.PhQWriteStruc( atoms.PosDir , true , System.Nsteps );
         }

	 // parameter output
	 FILE * Efile;  // energy file
	 Efile = fopen( "energy.out" , "w" );
	 WriteOutput MDSpec( Efile );

	 std::vector<Real> MDData;
	 MDData.resize( 4 );
	 
	 
	 //velocity init
	 Integrator.InitVelocities( atoms.VeloCar , System.dimension , System.Box , atoms.Mass );
	 atoms.Transform( atoms.VeloDir, atoms.VeloCar , System.RecBox );




	 System.Epot  =  MLFF.ComputeEnergy();
	 System.Ekin  =  ComputeKineticEnergy( atoms.VeloCar , atoms.Mass , atoms.Ekin );
	 System.Etot  =  System.Epot + System.Ekin ;
	 System.Tact  =  ComputeTemperature( System.Ekin , System.Natoms );

	 MDData = { System.Tact , System.Ekin , System.Epot , System.Etot };
	 //MDSpec.WriteVector( MDData );


	 // initialize analysis
	 PairDistFunction DistDistributions;
	 DistDistributions.DistCenter( System.Box[0][0] / Real( 6 ) , atoms.PosDir , 250 );
	 DistDistributions.PairDistInit( System.Box[0][0] , 250 );
	 AngularDist3d AngDist( 50 , 100 );


         std::vector<Real> COMVel;
         COMVel.resize( atoms.VeloDir[0].size() );


	 InsertParticle Insertion( 100 , System.temperature );


	 float t1 = omp_get_wtime();
	 for ( auto i = 0 ; i < System.Nsteps ; i++ ){
	 
	 	Integrator.ThermostatIntegrateMLFF( atoms.PosDir , atoms.ForceDir , atoms.ForceCar , atoms.VeloDir ,
	                          	        System.timeStep , System.Box , System.RecBox , MLFF , atoms.Mass );

	 	// switch to NVE ensemble
	 	if ( System.SwitchThermo == i ){
	 		 std::cout << "Switched to NVE ensemble at step " << i << std::endl;
	 		 System.Thermostyle = 0;
	 		 System.AndersenFrequ = Real( 0.0 );
	                 Integrator.InitAndersen( System.temperature , System.AndersenFrequ*System.timeStep , atoms.Mass );
	        }



	 	 //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ neighbor list update ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    	         atoms.Transform( atoms.VeloCar, atoms.VeloDir , System.Box );
	         atoms.Transform( atoms.PosCar, atoms.PosDir , System.Box );
	 	 MLFF.LinkList.UpdateNNVerletList( atoms.PosDir , atoms.PosCar , atoms.VeloCar ,
	 		                               System.timeStep , System.Cutoff , System.Box );
	 	// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

	 	// analysis and output
	 	if ( i%System.AnaFrequ == 0  and i > System.StartSample ){
	 	   //if nvt is used to equibrilate one can switch to nve after N steps

	 	    printf( "Evaluating sample Nr , %i \n" , i );
	            // compute pair dist
	 	    DistDistributions.ComputeDistToCenter( atoms.PosDir , System.Box );
	 	    AngDist.ComputeAngularCOMCentered( atoms.PosDir , DistDistributions.Center , System.Box );


	            atoms.Transform( atoms.ForceCar, atoms.ForceDir , System.Box );
	 	    //// compute md specs
                    MLFF.ComputeDescriptors( atoms.PosDir , System.Box );
	 	    System.Epot     =  MLFF.ComputeEnergy();
	            System.Ekin     =  ComputeKineticEnergy( atoms.VeloCar , atoms.Mass , atoms.Ekin );
	            System.Etot     =  System.Epot + System.Ekin ;
	            System.Tact     =  ComputeTemperature( System.Ekin , System.Natoms );
	 	    System.Pressure =  ComputePressure( atoms.PosDir , atoms.ForceCar ,
	 	    	                                System.Box , System.temperature ,
							MLFF.LinkList.NNList );

	 	    //// write md specs
	            MDData = { System.Tact , System.Ekin , System.Pressure , System.Epot, System.Etot };
	            MDSpec.WriteVector( MDData );
	 	    // write structure
	 	    if ( System.Output == 0 ){
	 	       StrucOut.AddStruc( atoms.PosDir );    //XDAT format
	 	    }
	 	    //Lammps
	 	    else if ( System.Output == 1 ){
	 	       StrucOut.LammpsWrite( atoms.PosDir , System.Box , atoms.Type ,
	 	                             System.Natoms );
	 	    }
	 	    if ( System.Output == 2 ){
	 	       StrucOut.PhQWriteStruc( atoms.PosDir , false , i );
	 	    }
	 	}
         }

	 float t2 = omp_get_wtime();

	 std::cout << "Execution time " << t2 - t1 << std::endl;

	 DistDistributions.WriteOutputCOM( "CenterDist.out" );
	 //DistDistributions.WriteOutputPD( "PairDist.out" );
	 AngDist.WriteOutput( "PolarDist.out" );

     return 0;	
 }


//  ++++++++++++++++++++++++++++++++++++++++++
//  ++++++++++++++++++++++++++++++++++++++++++++++++++
//
//  Harmonic potential with Higgs perturbation field
//
//  ++++++++++++++++++++++++++++++++++++++++++++++++++
//  +++++++++++++++++++++++++++++++++++++++++++


void MainStyles::SetUpHarmonicDoubleWellPerturbation( AtomType & atoms , SystemParams & System ,
                                                      std::string inputFile ){

	Units unit;

	Ran Rand( std::chrono::system_clock::to_time_t( std::chrono::system_clock::now() ) );
	
	FlagFinder Find( inputFile );

	System.dimension  =  int( Find.CheckFlag( "Ndim" ) );

	std::cout <<"##################################"<< std::endl;
	std::cout <<"########## PARAMETERS ############"<< std::endl;
	std::cout << "dimension " << System.dimension << std::endl;

	std::vector<Real> temp =  Find.CheckFlagVec( "NCells" );
	System.NParticles =  std::vector<int> ( temp.begin() , temp.end() );

	std::cout << "Repeating units " << System.NParticles[0] << "   "
		                        << System.NParticles[1] << "   "
                                        << System.NParticles[2] << std::endl;

	System.deltaX = Find.CheckFlagVec( "DeltaX" );
	std::cout << "Grid spacing [\\AA]" << std::setw( 8 ) << System.deltaX[0] << "   "
	   	                           << std::setw( 8 ) << System.deltaX[1] << "   "
                                           << std::setw( 8 ) << System.deltaX[2] << std::endl;


	System.Box = std::vector<std::vector<Real> > ( System.dimension , 
			          std::vector<Real> ( System.dimension , Real( 0.0 ) ) );
	System.RecBox = std::vector<std::vector<Real> > ( System.dimension , 
			          std::vector<Real> ( System.dimension , Real( 0.0 ) ) );

	System.timeStep = Find.CheckFlag( "TStep" ) * unit.timeConv;
	std::cout << "Time step in ps " << System.timeStep / unit.timeConv << std::endl;

	System.temperature = Find.CheckFlag( "Temp" );
	std::cout << "Temperature [K] " << System.temperature << std::endl;
	

        // determine cutoff radius in an Angstroem	
	System.Cutoff = Find.CheckFlag( "Cutoff" );
	std::cout << "Cutoff distance set to " << System.Cutoff << "[\AA]" << std::endl;


	std::cout << "Analysis frequency is set to" << std::endl;
	System.AnaFrequ =  int ( Find.CheckFlag( "AnaFrequ" ) );

	for ( auto i = 0 ; i < System.deltaX.size() ; i++ ){
		System.Box[ i ][ i ] = Real( System.NParticles[i] ) *
			               System.deltaX[i];
		System.RecBox[ i ][ i ] = Real( 1 ) / System.Box[ i ][ i ];
        }


	System.Natoms = std::accumulate( std::begin( System.NParticles ) ,
			                     std::end( System.NParticles ) , 1 ,
					     std::multiplies<Real>() );

	System.Nsteps = int( Find.CheckFlag( "NSTEPS" ) );

	std::cout << "Number of MD steps " << System.Nsteps << std::endl;


	atoms.PosDirData.resize( System.Natoms * System.dimension );
	atoms.InitializePointers( atoms.PosDirData , atoms.PosDir , System.dimension );

	atoms.PosCarData.resize( System.Natoms * System.dimension );
	atoms.InitializePointers( atoms.PosCarData , atoms.PosCar , System.dimension );

	atoms.VeloDirData.resize( System.Natoms * System.dimension );
	atoms.InitializePointers( atoms.VeloDirData , atoms.VeloDir , System.dimension );

	atoms.VeloCarData.resize( System.Natoms * System.dimension );
	atoms.InitializePointers( atoms.VeloCarData , atoms.VeloCar , System.dimension );

	atoms.ForceDirData.resize( System.Natoms * System.dimension );
	atoms.InitializePointers( atoms.ForceDirData , atoms.ForceDir , System.dimension );
	
	atoms.ForceCarData.resize( System.Natoms * System.dimension );
	atoms.InitializePointers( atoms.ForceCarData , atoms.ForceCar , System.dimension );

	atoms.Mass.resize( System.Natoms );
	atoms.Type.resize( System.Natoms );
	atoms.Spec.resize( System.Natoms );
	atoms.Epot.resize( System.Natoms );
	atoms.Ekin.resize( System.Natoms );

	// reading thermostat stuff
	System.Thermostyle  =  int( Find.CheckFlag( "ThermST" ) );
	if ( System.Thermostyle == 0 ){
		std::cout << "Andersen thermostat is used" << std::endl;
		System.AndersenFrequ  =  Find.CheckFlag( "COLFA" ) / unit.timeConv ;
		std::cout << "Andersen collision frequ " << std::setw( 8 ) <<
		System.AndersenFrequ * System.timeStep <<std::endl;
        }
	else if ( System.Thermostyle == 1 ){
		std::cout << "Nose-Hoover thermostat is used" << std::endl;
		System.NoseMass  =  Find.CheckFlag( "NoseM" ) * unit.timeConv ;
		std::cout << "Nose-Hoover mass " << std::setw( 8 ) <<
		              System.NoseMass << std::endl;
        } 
	else if ( System.Thermostyle == 2 ){
		std::cout << "Langevin thermostat is used" << std::endl;
		System.LangevinFriction  =  Find.CheckFlag( "LangFric" ) * unit.timeConv ;
		std::cout << "Langevin friction coefficient " << std::setw( 8 ) <<
		              System.LangevinFriction << std::endl;
        }
	else if ( System.Thermostyle == 3 ){
		std::cout << "Bussi-Donadio-Parrinello thermostat is used" << std::endl;
		System.BDPThermoTau  =  Find.CheckFlag( "BDPTau" ) * unit.timeConv ;
		std::cout << "BDP time scale tau " << std::setw( 8 ) <<
		              System.BDPThermoTau / unit.timeConv << std::endl;
	}




        // delta vector for phonon dispersion in
	// harmonic approximation
	// supply in cartesian coordinates	
	System.DeltaPhonon = Find.CheckFlagVec( "DeltaPhon" );
	if ( System.DeltaPhonon.size() == 0 ){
		System.PhononOnOff = false;
        }
	else{
		System.PhononOnOff = true;
		System.DeltaPhonon = MatVec( System.RecBox , System.DeltaPhonon );
	}


	System.Output  =  0;
	System.Output  =  int( Find.CheckFlag( "Output" ) );
	if ( System.Output < 0 ){
		System.Output = 0;
        }


	Real M = Find.CheckFlag( "mass" );
	std::cout << "Particle masses [Da] " << std::setw( 8 ) <<
		         M <<std::endl;


	System.StartSample = int( Find.CheckFlag( "EquiStep" )  );

	if ( System.StartSample > 0 ){
		std::cout <<  "Starting to take samples at MD step "
			      <<  System.StartSample << std::endl;
        }else{
		std::cout << "No equibrilation time set; start sampling at first step" << std::endl;
		System.StartSample = -1;
        }

	System.SwitchThermo = int( Find.CheckFlag( "SwitchNVE" ) );
	if ( System.SwitchThermo > 0 ){
	      std::cout << "Switching to NVE ensemble after " << 
	      System.SwitchThermo << " steps activated" << std::endl;
	}



	System.NNHarmonicParams  = Find.CheckFlagVec( "NNHarmonic" );
	if ( System.NNHarmonicParams.size() == 0 ){
		std::cout << "No parameters for the harmonic potential supplied" << std::endl;
		std::cout << "Nobody coded default parameters; so I set everything" << std::endl;
		std::cout << "to zero. This will be an ideal gas" << std::endl;
		System.NNHarmonicParams.push_back( 0.0 );
		System.NNHarmonicParams.push_back( 0.0 );
	}
	else{
		std::cout << "Harmonic parameters set to " << std::endl;
		for ( auto i = 0 ; i < System.NNHarmonicParams.size() ; i++ ){
			std::cout << System.NNHarmonicParams[i] << std::endl;
		}
	}
	
	
	System.ONHarmonicParams  = Find.CheckFlagVec( "ONHarmonic" );
	if ( System.NNHarmonicParams.size() == 0 ){
		std::cout << "No parameters for the harmonic potential supplied" << std::endl;
		std::cout << "Nobody coded default parameters; so I set everything" << std::endl;
		std::cout << "to zero. This will be an ideal gas" << std::endl;
		System.ONHarmonicParams.push_back( 0.0 );
	}
	else{
		std::cout << "ON SITE Harmonic parameters set to " << std::endl;
		for ( auto i = 0 ; i < System.ONHarmonicParams.size() ; i++ ){
			std::cout << System.ONHarmonicParams[i] << std::endl;
		}
	}



	System.NNDoubleWellParams  =  Find.CheckFlagVec( "NNDoubleWell" );
	if ( System.NNDoubleWellParams.size() == 0 ){
		std::cout << "No parameters for the double-well potential supplied" << std::endl;
		std::cout << "Nobody coded default parameters; so I set everything" << std::endl;
		std::cout << "to zero. This will be an ideal gas" << std::endl;
                System.NNDoubleWellParams.push_back( 0.0 );
                System.NNDoubleWellParams.push_back( 0.0 );
                System.NNDoubleWellParams.push_back( 0.0 );
	}
	else{
		std::cout << "Double well parameters set to " << std::endl;
		for ( auto i = 0 ; i < System.NNDoubleWellParams.size() ; i++ ){
			std::cout << System.NNDoubleWellParams[i] << std::endl;
		}
	}

	std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;




	System.RestartFile =  Find.CheckStringFlag( "restart_file" );
	unsigned int len   =  System.RestartFile.length();

	if ( !System.RestartFile.substr( 0, len ).compare( "NOTFOUND" ) ){
	    unsigned int col = 0;
	    for ( auto i = 0 ; i < System.NParticles[0] ; i++ ){
	      	Real x = Real( i ) * System.deltaX[0] + System.deltaX[0] * Real( 0.5 );
	       	for ( auto j = 0 ; j < System.NParticles[1] ; j++ ){
	       	    Real y = Real( j ) * System.deltaX[1] + System.deltaX[1] * Real( 0.5 );
	       		for ( auto k = 0 ; k < System.NParticles[2] ; k++ ){
	       			Real z = Real( k ) * System.deltaX[2] + System.deltaX[2] * Real( 0.5 );
	       			atoms.Type[ col ] = 1;
		   		    if ( !System.PhononOnOff ){
	       			            atoms.PosCar[ col ][ 0 ] = x ;
	       			            atoms.PosCar[ col ][ 1 ] = y ;
	       			            atoms.PosCar[ col ][ 2 ] = z ;
		   		    }
		   		    else{
					    atoms.PosCar[ col ][ 0 ] = x ;
	       			            atoms.PosCar[ col ][ 1 ] = y ;
	       			            atoms.PosCar[ col ][ 2 ] = z ;
		   		   }
	       			   atoms.PosDir[ col ] = MatVec( System.RecBox , atoms.PosCar[ col ] );
	       			   atoms.Mass[col] = M;
	       			   col++;
	       	    }
	       	}
	    }
        }else{
		  atoms.ReadRestartStructure( System.RestartFile );
		  std::cout << "Old structure read " << std::endl;
		  std::cout << "Check that your other parameters are in" << std::endl;
		  std::cout << "agreement with the supplied structure" << std::endl;
		  for ( auto i = 0 ; i < atoms.PosDir.size() ; i++ ){
			  std::cout << atoms.PosDir[i][0] << "   "
		  	            << atoms.PosDir[i][1] << "    "
	   	     	    	    << atoms.PosDir[i][2] << std::endl;
	          atoms.Mass[ i ] = M;
		  }
	}


	System.ParticleInsert = ( int ) ( Find.CheckFlag( "ParticleInsert" )  );
	if ( System.ParticleInsert > 0 ){
		std::cout << "New particle will be inserted every "
			       << System.ParticleInsert << " steps" << std::endl;

	}else{
		std::cout << "Particle insertion switched off" << std::endl;
	}

	
	// read heat transport variables 
	System.HeatVars.HeatGradient = ( int ) ( Find.CheckFlag( "HeatGrad" )  );
	if ( System.HeatVars.HeatGradient > 0 ){
		std::cout << "Heat gradient was switched on" << std::endl;
	        System.HeatVars.temperatureA =  Find.CheckFlag( "TempA" );
	        System.HeatVars.temperatureB =  Find.CheckFlag( "TempB" );
		if ( System.HeatVars.temperatureA < 0 ){
			throw ( System.HeatVars.temperatureA );
		}
		std::cout << "Temperature A for heat transport set to" 
			  <<  std::setw( 8 ) << System.HeatVars.temperatureA << std::endl;
		if ( System.HeatVars.temperatureB < 0 ){
			throw ( System.HeatVars.temperatureB );
		}
		std::cout << "Temperature B for heat transport set to" 
			  <<  std::setw( 8 ) << System.HeatVars.temperatureB << std::endl;

	        System.HeatVars.Direction    =  ( int ) ( Find.CheckFlag( "HeatDir" ) );
		if ( System.HeatVars.Direction < 1 or System.HeatVars.Direction > 3 ){
			std::cout << "Direction for temperature gradient not given" << std::endl;
			std::cout << "Supply a number {1,2,3}" << std::endl;
			std::cout << "x ---> 1" << std::endl;
			std::cout << "y ---> 2" << std::endl;
			std::cout << "z ---> 3" << std::endl;
		}
		else{
			std::cout << "HeatGardient direction set to " 
				  << System.HeatVars.Direction << std::endl;
		}

		System.HeatVars.ActiveMeasure = ( int ) ( Find.CheckFlag( "ActMeasure" ) );
		if ( System.HeatVars.ActiveMeasure < 0 ){
			std::cout << "Active heat gradient measure switched off" << std::endl;
		}
		else{
			std::cout << "Active heat gradient measure switched on" << std::endl;
		}

		System.HeatVars.EquiSteps  =  ( int ) ( Find.CheckFlag( "ActMEqui" ) );
	}

	std::cout <<"##################################"<< std::endl;
}


 

int MainStyles::HarmonicDoubleWellPerturbation( const char * fname ){

	AtomType atoms;       // contains information about atoms
	SystemParams System;  // information about the simulation conditions

	// reading input parameters
	SetUpHarmonicDoubleWellPerturbation( atoms , System , fname );



	//
	// if grad descent 
	// init grad descent here 
	// 
	//
	//



	// set up force class
	ComputeForces Force( atoms.PosDir , ( unsigned int )( 6 ) , System.Box );

	// set up harmonic force  coupling strength equilibrium dist cart equilibrium dist dir
	Force.SetUpHarmonicPart( System.NNHarmonicParams[1] ,
			         System.NNHarmonicParams[0] );

	// set up the double well potential coupling parameter of second order
	//                                  coupling parameter of fourth order
	//                                  equilibrium distance in direct coordinates
	//
	Force.SetUpDoubleWellPart( System.NNDoubleWellParams[0] ,
			           System.NNDoubleWellParams[1] ,
				   System.NNDoubleWellParams[2] ,
			           System.NNDoubleWellParams[3] );


	Force.SetUpOnSitePotential( System.ONHarmonicParams[0] , 
			    int( System.ONHarmonicParams[1] ) , atoms.PosDir );


        ThermostatsVeloVerlet Integrator;
	// set up thermostat integrator velocity verlet
	if ( System.HeatVars.HeatGradient < 1 ){
	    if ( System.Thermostyle == 0 ){
	            Integrator.InitAndersen( System.temperature ,
				System.AndersenFrequ*System.timeStep , atoms.Mass );
	    }
	    else if ( System.Thermostyle == 1 ){
	    	Integrator.InitNoseHoover( System.NoseMass , 
	          	   atoms.PosDir.size() * atoms.PosDir[0].size() + 1 ,
			   System.temperature , atoms.Mass );
	    }
	    else if ( System.Thermostyle == 2 ){
	    	Integrator.InitLangevin( System.LangevinFriction ,
				         System.temperature , 
	    			         System.timeStep , System.dimension , atoms.Mass );
	    }
	    else if ( System.Thermostyle == 3 ){
		    Integrator.InitBDP( System.BDPThermoTau , System.temperature ,
				        System.timeStep , 
					atoms.PosDir.size() * atoms.PosDir[0].size() ,
					// total degrees of freedom
					atoms.Mass );
	    }
	}
	// do the heat gradient initialization
	else{
	    if ( System.Thermostyle == 0 ){
		System.HeatVars.nu = System.AndersenFrequ * System.timeStep;
		Integrator.InitHeatGradientThermoAndersen( System.HeatVars , 
				System.NParticles , atoms.Mass );
	    }
	    else if ( System.Thermostyle == 2 ){
	    	Integrator.InitHeatGradientThermoLangevin( 
				         System.HeatVars , System.LangevinFriction ,
	    			         System.timeStep , System.dimension , 
					 System.NParticles , atoms.Mass );
	    }
	    else if ( System.Thermostyle == 3 ){
		    Integrator.InitHeatGradientThermoBDP( System.HeatVars , System.BDPThermoTau ,
				System.timeStep , atoms.Mass , System.NParticles , 
				System.NParticles.size() );
	    }
	}

	// force initialization function call has to match
	// the one used for force computation in andersen routine
	//atoms.ForceDir =  Force.SupplyHarmonicOnly( atoms.PosDir );
	//Force.SupplyMorseOnly( atoms.PosDir , atoms.ForceDir , System.Box , System.RecBox );
	//atoms.ForceDir = Force.SupplyHarmonicAndOnSite( atoms.PosDir );



	//phonon finite differences is computed here
	if ( System.PhononOnOff ){
		FDPhonons Phonons( System.Box , System.RecBox , System.NParticles , System.DeltaPhonon );
		Phonons.main( atoms.PosDir , Force , System.Box , atoms.Mass );
		return 0;
	}



	std::vector<std::string> types;
	types = { "H" };
	std::vector<int> Nats ( 1 , System.Natoms );

	// STRUCTURE output
	FILE * xdatcar;
	if ( System.Output == 0 ){
		xdatcar = fopen( "XDATCAR" , "w" );
	}else if ( System.Output == 1 ){
		xdatcar = fopen( "traj.lammpstrj" , "w" );
        }
	else if ( System.Output == 2 ){
		xdatcar = fopen( "md.out" , "w" );
	}

	//FILE * ForceCar;
	//ForceCar = fopen( "Forces" , "w" );


        WriteOutput StrucOut( xdatcar );
        //WriteOutput ForceOut( ForceCar );
	// xdat init
	//StrucOut.XdatInit( System.Box , types , Nats );
	//StrucOut.AddStruc( atoms.PosDir );

	if ( System.Output == 0 ){
		StrucOut.XdatInit( System.Box , types , Nats );
	        StrucOut.AddStruc( atoms.PosDir );
        }else if ( System.Output == 1 ){
            	StrucOut.LammpsWrite( atoms.PosDir , System.Box , atoms.Type ,
            	               System.Natoms );
        }
            else if ( System.Output == 2 ){
            	StrucOut.PhQWriteStruc( atoms.PosDir , true , System.Nsteps );
        }

	// parameter output
	FILE * Efile;  // energy file
	Efile = fopen( "energy.out" , "w" );
	WriteOutput MDSpec( Efile );

	std::vector<Real> MDData;
	MDData.resize( 4 );


	//velocity init
	Integrator.InitVelocities( atoms.VeloCar , System.dimension , System.Box , atoms.Mass );
	atoms.Transform( atoms.VeloDir, atoms.VeloCar , System.RecBox );



        //System.Epot  =  Force.PotentialEnergyHarmonic( atoms.PosDir , System.Box );
        //System.Epot  =  Force.PotentialEnergyMorse( atoms.PosDir , System.Box );
        //System.Epot  =  Force.PotentialEnergyHarmonicAndOnSite( atoms.PosDir , System.Box );
	System.Epot    =  Force.PotentialEnergyDoubleWell( atoms.PosDir , System.Box );
	System.Ekin  =  ComputeKineticEnergy( atoms.VeloCar , atoms.Mass , atoms.Ekin );
	System.Etot  =  System.Epot + System.Ekin ;
	System.Tact  =  ComputeTemperature( System.Ekin , System.Natoms );

	MDData = { System.Tact , System.Ekin , System.Epot , System.Etot };
	MDSpec.WriteVector( MDData );


	// initialize analysis
	PairDistFunction DistDistributions;
	DistDistributions.DistCenter( System.Box[0][0] / Real( 6 ) , atoms.PosDir , 250 );
	DistDistributions.PairDistInit( System.Box[0][0] , 10000 );
	AngularDist3d AngDist( 50 , 100 );
	

	AnalyzePotential PotentialAna( atoms.PosDir );


        std::vector<Real> COMVel;
        COMVel.resize( atoms.VeloDir[0].size() );




	// this part can be used to test the potential energy computation

//	{
//
//	   int Steps    =  100;
//	   Real dstep   =  atoms.PosDir[1][2] - atoms.PosDir[0][2];
//	   Real EquiPos =  atoms.PosDir[0][2] * System.Box[0][0];
//           dstep        =  dstep  /  Real( Steps );
//
//	   atoms.PosDir[0][2]   =  atoms.PosDir[0][2] - dstep * Real( Steps ) / Real( 2 );
//
//
//	   FILE * Efile;  // energy file
//	   Efile = fopen( "test.out" , "w" );
//
//
//	   for ( auto i = 0 ; i < Steps ; i++ ){
//	   	//System.Epot    =  Force.PotentialEnergyDoubleWell( atoms.PosDir , System.Box );
//	   	//Force.SupplyDoubleWellOnly( atoms.PosDir , System.Box , System.RecBox , atoms.ForceDir );
//		System.Epot     = Force.SupplyOnSiteEnergyOnly( atoms.PosDir ,
//		  	                                        System.Box );
//		Force.SupplyOnSiteForceOnly( atoms.PosDir , System.Box , 
//				             System.RecBox , atoms.ForceDir );
//		
//		fprintf( Efile , "   %f" , atoms.PosDir[0][2] * System.Box[0][0] - EquiPos );
//		fprintf( Efile , "   %f" , System.Epot );
//		fprintf( Efile , "   %f" , atoms.ForceDir[0][0] * System.Box[0][0] );
//		fprintf( Efile , "   %f" , atoms.ForceDir[0][1] * System.Box[0][0] );
//		fprintf( Efile , "   %f" , atoms.ForceDir[0][2] * System.Box[0][0] );
//		fprintf( Efile , "\n" );
//                fflush( Efile );
//	   	atoms.PosDir[0][2]  =  atoms.PosDir[0][2] + dstep;
//	   }
//
//	   throw;
//	}








	for ( auto i = 0 ; i < System.Nsteps ; i++ ){


		Integrator.ThermostatIntegrateHarmonicDoubleWellPerturbation(
				atoms.PosDir , atoms.ForceDir , atoms.VeloDir ,
	                        System.timeStep , System.Box , 
				System.RecBox , Force , atoms.Mass );

		if ( System.SwitchThermo == i ){
			std::cout << "Switched to NVE ensemble at step " << i << std::endl;
			System.Thermostyle = 0;
			System.AndersenFrequ = Real( 0.0 );
          	        Integrator.InitAndersen( System.temperature , 
					System.AndersenFrequ*System.timeStep , atoms.Mass );
	        }

		// analysis and output
		if ( i%System.AnaFrequ == 0  and i > System.StartSample ){
		   //if nvt is used to equibrilate one can switch to nve after N steps

		   printf( "Evaluating sample Nr , %i \n" , i );
	           // compute pair dist
		   //DistDistributions.ComputeDistToCenter( atoms.PosDir , System.Box );
		   //DistDistributions.ComputePairDistribution( atoms.PosDir , System.Box );
		   //AngDist.ComputeAngularCOMCentered( atoms.PosDir , DistDistributions.Center , System.Box );

	           atoms.Transform( atoms.VeloCar, atoms.VeloDir , System.Box );

		   // compute md specs
	 	   //System.Epot   =   Force.PotentialEnergyDoubleWell( atoms.PosDir , System.Box );
		   //System.Epot   =   Force.PotentialEnergyHarmonic( atoms.PosDir , System.Box );
		   //System.Epot     = Force.SupplyOnSiteEnergyOnly( atoms.PosDir ,
		   //		                                   System.Box );
		   System.Epot   =   Force.PotentialEnergyHarmonicAndOnSite( atoms.PosDir ,
				                                             System.Box );
	           System.Ekin   =   ComputeKineticEnergy( atoms.VeloCar , atoms.Mass , atoms.Ekin );
	           System.Etot   =   System.Epot + System.Ekin ;
	           System.Tact   =   ComputeTemperature( System.Ekin , System.Natoms );
		   
		   // compute layer temperatures and write them to thermostats
		   // debug file
		   

		   // write md specs
	           //MDData = {System.Tact , System.Ekin , System.Epot, System.Etot , COMVel[0],COMVel[1],COMVel[2] };
	           MDData = {System.Tact , System.Ekin , System.Epot, System.Etot };
	           MDSpec.WriteVector( MDData );
		   // write structure
		   if ( System.Output == 0 ){
		           StrucOut.AddStruc( atoms.PosDir );    //XDAT format
		   }
		   //Lammps
		   else if ( System.Output == 1 ){
		           StrucOut.LammpsWrite( atoms.PosDir , System.Box , atoms.Type ,
		                         System.Natoms );
		   }
		   if ( System.Output == 2 ){
		           StrucOut.PhQWriteStruc( atoms.PosDir , false , i );
		   }
		}
        }
	//DistDistributions.WriteOutputCOM( "CenterDist.out" );
	DistDistributions.WriteOutputPD( "PairDist.out" );
	//AngDist.WriteOutput( "PolarDist.out" );
	//PotentialAna.PotentialAnaFinalize();

	
    return 0;	
}

