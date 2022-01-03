

#ifndef _ANALYZE
#include "analyze.h"
#endif

 unsigned int computeBin( Real x , Real dx , unsigned int Max ){

	 unsigned int Bin = floor( x / dx );
	 if ( Bin > Max - 1 ){ Bin = Max - 1; }
	 if ( Bin < 0 ){ Bin =  0;}
	 return Bin;
 }



 Real compute_spherical_theta( std::vector<Real> x ){

	 Real theta ;

	 if ( x.size() != 3 ){
		 printf( "~~~~~~~~   ERROR ~~~~~~~~~~~~\n" );
		 printf( "Error in function compute_spherical_theta\n" );
		 printf( "Size of vector does not fit\n" );
		 printf( "Only 3-vectors are allowed\n" );
		 printf( "Returning zero\n" );
		 printf( "~~~~~~~~   ERROR ~~~~~~~~~~~~\n" );
		 theta = Real( 0 );
		 return theta;
	 }

	 Real norm = VecNorm( x );


         if ( std::abs( x[ 2 ] ) > Real( 1e-12 ) ){
           theta  =  atan( sqrt( ( x[ 0 ] / norm )*( x[ 0 ] / norm )
			       + ( x[ 1 ] / norm )*( x[ 1 ] / norm ) )
                                 / ( x[ 2 ] / norm ) );
         }
         else{
            theta  =  M_PI / Real( 2 );
	 }

         if ( x[ 2 ] < Real( 0 ) ){
            theta  =  M_PI + theta;
	 }

         if ( theta < Real( 0 ) ){
            theta  =  std::abs( theta );
         }

         if ( theta > M_PI ){
           theta  =  theta - M_PI;
	 }

	 return theta;


 }

 Real compute_spherical_phi( std::vector<Real> x ){

	 Real phi;

/*	 Print( "" );
	 Print( "" );
	 Print( "" );
	 Print( "" );
	 */
	 if ( x.size() != 3 ){
		 printf( "~~~~~~~~   ERROR ~~~~~~~~~~~~\n" );
		 printf( "Error in function compute_spherical_phi\n" );
		 printf( "Size of vector does not fit\n" );
		 printf( "Only 3-vectors are allowed\n" );
		 printf( "Returning zero\n" );
		 printf( "~~~~~~~~   ERROR ~~~~~~~~~~~~\n" );
		 phi = Real( 0 );
		 return phi;
	 }

	 //Print( x );

	 Real norm = VecNorm( x );
	 phi = atan( ( x[1] / norm ) / ( x[0] / norm ) );

	  if ( std::abs( x[ 0 ] ) <  Real( 1e-7 ) && std::abs( x[1] ) < Real( 1e-7 ) ){
             phi  =  Real( 0 );
	  }
	  else if ( x[ 0 ] > Real( 0 ) && x[ 1 ] < Real( 0 ) ){
            phi  =  Real( 2 ) * M_PI + phi;
          }
	  else if ( std::abs( x[0] ) < Real( 1e-7 ) && x[1] > Real( 0 ) ){
            phi  =  M_PI * Real( 0.5 );
	  }
	  else if ( std::abs( x[0] ) < Real( 1e-7 ) && x[ 1 ] < Real( 0 ) ){
            phi  =  Real( 3 ) * M_PI / Real( 2 );
	  }
          else if ( x[ 0 ] < Real( 0 ) ){
            phi  =  phi + M_PI;
          }

	  if ( phi < Real( 0 ) ){
		  phi = M_PI * Real( 2 ) + phi;
	  }
	  else if( phi > Real( 2 ) * M_PI ){
		  phi = phi - Real( 2 ) * M_PI;
	  }

	  return phi;

 }


 void SphericalSinCorrection( std::vector<std::vector<Real> > & data , int dim ){

	 Real theta;
	 if ( dim == 0 ){
	   Real dx  = M_PI / Real( data.size() );
	   for ( auto i = 0 ; i < data.size() ; i++ ){
		   theta = Real( i ) * dx + dx * Real( 0.5 );
		   for ( auto j = 0 ; j < data[i].size(); j++ ){
			   data[i][j]  = data[i][j] / sin( theta );
	           }
           }
         }
	 else{
	   Real dx  = M_PI / Real( data[0].size() );
	   for ( auto i = 0 ; i < data.size() ; i++ ){
		   for ( auto j = 0 ; j < data[i].size(); j++ ){
		           theta = Real( j ) * dx + dx * Real( 0.5 );
			   data[i][j]  = data[i][j] / sin( theta );
	           }
           }

         }


 }



 Real ComputeKineticEnergy( std::vector<VecNSlice<Real> > Velo , std::vector<Real> Mass , std::vector<Real>& Ekin ){

	Real Ek = Real( 0 );
//#ifdef USE_OPENMP
//#pragma omp parallel for reduction(+:Ek)
//#endif
	for ( auto i = 0 ; i < Velo.size() ; i++ ){
		Ekin[i] = Real( 0.5 ) * DotProduct( Velo[i] , Velo[i] ) * Mass[i];
		Ek = Ek + Ekin[ i ] ;
	}
	return Ek;
 }


 Real ComputeTemperature( Real Ek , int Nf ){
	 Units units;
	 Real T =  Real( 2 ) * Ek / Real( 3 ) / Real( Nf ) / units.kB;
	 return T;
 }



 // compute pressure in NVT ensemble
 Real ComputePressure( std::vector<VecNSlice<Real> > Pos , std::vector<VecNSlice<Real> > force ,
		               const std::vector<std::vector<Real> > Box , const Real T ,
					   const std::vector<std::vector<int> > List ){

	 Real pressure = Real( 0 );
	 Units units;
	 Real volume = DotProduct( Box[ 2 ] , CrossProduct3D( Box[ 0 ] , Box[ 1 ] ) );

	 Real rho = Real( Pos.size() ) / volume;
	 Real t1 = rho * T * units.kB;
	 Real Mean = Real( 0 );
//#ifdef USE_OPENMP
//#pragma omp parallel for reduction(+:Mean)
//#endif
	 for ( auto i = 0 ; i < Pos.size() ; i++ ){
		 for ( auto j = 0 ; j < List[i].size() ; j++ ){
			 std::vector<Real> delta =   MatVec( Box , get_nearest_image( Pos[ i ] , Pos[ List[ i ][ j ] ] ) - 
				                         Pos[ i ] );
			 Mean = Mean + DotProduct( force[i] , delta ); 
		 }
	 }
	 pressure = t1 + Real( 0.5 ) / ( Real( Pos[0].size() ) * volume ) * Mean / Real( Pos.size() ); 
	 return pressure;
 }






 /*
  *
  *Pair distribution function
  *
  */


 void PairDistFunction::DistCenter( Real xmax , std::vector<VecNSlice<Real> > Cent , int Nmax ){

	 distribution.resize( Nmax );
	 xaxis.resize( Nmax );
	 dx = xmax / Real( Nmax );
	 Nbins = Nmax;
	 for ( auto i = 0 ; i < Nmax ;i++ ){
		 xaxis[i] = dx * Real( i ) + dx / Real( 2 );
         }
	 Center.resize( Cent.size() );
	 for ( auto i = 0 ; i < Cent.size() ; i++ ){
		 Center[i].resize( Cent[i].size() );
		 for ( auto j = 0 ; j < Cent[i].size(); j++ ){
			 Center[i][j] = Cent[i][j];
		 }
	 }
 }

void PairDistFunction::ComputeDistToCenter( std::vector<VecNSlice<Real> > Pos ,
		                                 std::vector<std::vector<Real>> lattice ){

	 for ( auto i = 0 ; i < Pos.size() ; i++ ){
		 std::vector<Real> TempPos = get_nearest_image( Center[i] , Pos[i] );
		 std::vector<Real> delta = MatVec( lattice , SubVecs( TempPos , Center[i] ) );
		 Real norm = VecNorm( delta );
	 	 unsigned int bin = computeBin( norm , dx , Nbins );
		 distribution[bin] = distribution[bin] + Real( 1 );
     }
 }



 void PairDistFunction::PairDistInit( Real xmax , int Nmax ){

	 PairDistribution.resize( Nmax );
	 PairAxis.resize( Nmax );
	 PairDx = xmax / Real( Nmax );
	 PairBins = Nmax;
	 for ( auto i = 0 ; i < Nmax ;i++ ){
		 PairAxis[i] = PairDx * Real( i ) + PairDx / Real( 2 );
     }
 }

 void PairDistFunction::ComputePairDistribution( std::vector<VecNSlice<Real> > Pos , std::vector<std::vector<Real> > Box ){

//#pragma omp parallel for reduction(vec_Real_plus : PairDistribution )
	 for ( auto i = 0 ; i < Pos.size() - 1 ; i++ ){
		 for ( auto j = i + 1 ; j < Pos.size() ; j++ ){
			 std::vector<Real> r = get_nearest_image( Pos[ i ] , Pos[j] ) - Pos[i];
			 Real norm = VecNorm( MatVec( Box , r ) );
			 unsigned int Bin = computeBin( norm , PairDx , PairBins );
			 PairDistribution[ Bin ] =  PairDistribution[ Bin ] + Real( 1 ); 
		 }
	 }
 }



void PairDistFunction::WriteOutputCOM( std::string fname ){


	Real norm = D1_simps( distribution );
	norm = Real( 1 ) / ( norm * dx ) ;
        distribution = ScaTVec( norm , distribution );

	FILE * outfile;  // energy file
        outfile = fopen( fname.c_str() , "w" );
	for ( auto i = 0 ; i < distribution.size() ; i++ ){
		fprintf( outfile , "%f    %f \n" , xaxis[i] , distribution[i] );
	}

	fclose( outfile );
 }


void PairDistFunction::WriteOutputPD( std::string fname ){


	Real norm = D1_simps( PairDistribution );
	norm = Real( 1 ) / ( norm * PairDx ) ;
    PairDistribution = ScaTVec( norm , PairDistribution );

	FILE * outfile;  // energy file
        outfile = fopen( fname.c_str() , "w" );
	for ( auto i = 0 ; i < PairDistribution.size() ; i++ ){
		fprintf( outfile , "%f    %f \n" , PairAxis[i] , PairDistribution[i] );
	}

	fclose( outfile );
 }



 /*
  *
  * Angluar distrinution in spherical coords
  *
  */


 AngularDist3d::AngularDist3d( unsigned int N1 , unsigned int N2 ){

	 Ntheta  =  N1;
	 Nphi    =  N2;
	 dtheta  =  M_PI / Real( N1 );
	 dphi    =  Real( 2 ) * M_PI / Real( N2 );
	 distribution.resize( N1 );
	 for ( auto i = 0 ; i < N1 ; i++ ){
	     distribution[i].resize( N2 );
	 }
  }

 void AngularDist3d::ComputeAngularCOMCentered( std::vector<VecNSlice<Real>> Pos ,
		                                std::vector<std::vector<Real>> COM ,
					        std::vector<std::vector<Real>> lattice ){

	 for ( auto i = 0 ; i < Pos.size() ; i++ ){
		 std::vector<Real> TempPos = get_nearest_image( COM[i] , Pos[i] );
		 std::vector<Real> delta = MatVec( lattice , SubVecs( TempPos , COM[i] ) );

		 Real theta =  compute_spherical_theta( delta );
		 Real phi   =  compute_spherical_phi( delta );

		 unsigned int thetaBin = computeBin( theta , dtheta , Ntheta );
		 unsigned int phiBin = computeBin( phi , dphi , Nphi );
		 distribution[ thetaBin ][ phiBin ] =
			    distribution[ thetaBin ][ phiBin ] + Real( 1 );
	 }
 }

 void AngularDist3d::WriteOutput( std::string fname ){

	SphericalSinCorrection( distribution );
	Real norm = D2_simps( distribution );
	norm = Real( 1 ) / ( norm * dtheta * dphi ) ;
        distribution = ScaTMat( norm , distribution );

	FILE * outfile;  // energy file
        outfile = fopen( fname.c_str() , "w" );

	for ( auto i = 0 ; i < Ntheta ; i++ ){
	     for ( auto j = 0 ; j < Nphi ; j++ ){
		     fprintf( outfile , "    %f" , distribution[i][j] );
	     }
	     fprintf( outfile , "\n" );
        }

	fclose( outfile );
 }



AnalyzePotential::AnalyzePotential( const std::vector<VecNSlice<Real> >& Pos ){


	unsigned int dim = 3;
        EquiPosArray.resize( dim * Pos.size() );
        unsigned int col = 0;
        for ( size_t i = 0 ; i < Pos.size() ; i++ ){
                for ( size_t j = 0 ; j < Pos[i].size() ; j++ ){
                        EquiPosArray[col] = Pos[i][j];
                        col++;
                }
        }

        AssignVector( EquiPosArray , EquiPos , dim );

	CenterArray.resize( dim );
	for ( size_t i = 0 ; i < dim ; i++ ){
		CenterArray[i] = EquiPos[0][i];
	}
	
	VecNSlice<Real> A( CenterArray.data() , dim );
	Center.push_back( A );

	DeltaArray.resize( dim );
	VecNSlice<Real> B ( DeltaArray.data() , dim );
	Delta.push_back( B );

	MakeAxisDirections();
	MakeFaceDirections();
	MakeRoomDirections();


        // first entry belongs to axis displacement
        // second entry belongs to face displacements
        Histograms.resize( 3 );
        Counter.resize( 3 );
        ForceHistogram.resize( 3 );
        for ( size_t i = 0 ; i < 3 ; i++ ){
                Histograms[i].resize( Nbins );
                ForceHistogram[i].resize( Nbins );
                Counter[i].resize( Nbins );
        }
}


void AnalyzePotential::MakeAxisDirections( void ){
	AxisArray.resize( 3*6 );
	AxisArray[0] = Real(1);
	AxisArray[4] = Real(1);
	AxisArray[7] = Real(1);
	AxisArray[9] = Real(-1);
	AxisArray[13] = Real(-1);
	AxisArray[17] = Real(-1);
	for ( size_t i = 0 ; i < AxisArray.size() ; i=i+3 ){
		VecNSlice<Real> A ( AxisArray.data()+i , 3 );
		Axis.push_back( A );
	}
}


void AnalyzePotential::MakeFaceDirections( void ){
	FaceArray.resize( 3*12 );
	FaceArray[0] = 1;
	FaceArray[1] = 1;
	FaceArray[3] = -1;
	FaceArray[4] = 1;
	FaceArray[6] = 1;
	FaceArray[7] = -1;
	FaceArray[9] = -1;
	FaceArray[10] = -1;

	FaceArray[12] = 1;
	FaceArray[14] = 1;
	FaceArray[15] = -1;
	FaceArray[17] = 1;
	FaceArray[18] = 1;
	FaceArray[20] = -1;
	FaceArray[21] = -1;
	FaceArray[23] = -1;

	FaceArray[25] = 1;
	FaceArray[26] = 1;
	FaceArray[28] = -1;
	FaceArray[29] = 1;
	FaceArray[31] = 1;
	FaceArray[32] = -1;
	FaceArray[34] = 1;
	FaceArray[35] = -1;




	unsigned int col = 0;
	for ( size_t i = 0 ; i < FaceArray.size() ; i=i+3 ){
		VecNSlice<Real> A ( FaceArray.data()+i , 3 );
		Face.push_back( A );
		Real Norm = sqrt( DotProduct( Face[col] , Face[col] ) );
		Face[col] = Face[col] / Norm;
		col++;
	}
}


void AnalyzePotential::MakeRoomDirections( void ){

	RoomArray.resize( 3*8 );
	RoomArray[0] = 1;
	RoomArray[1] = 1;
	RoomArray[2] = 1;

	RoomArray[3] = -1;
	RoomArray[4] = 1;
	RoomArray[5] = 1;
	
	RoomArray[6] = 1;
	RoomArray[7] = -1;
	RoomArray[8] = 1;
	
	RoomArray[9]  = 1;
	RoomArray[10] = 1;
	RoomArray[11] = -1;

	RoomArray[12] =  -1;
	RoomArray[13] =  -1;
	RoomArray[14] =   1;

	RoomArray[15] =  -1;
	RoomArray[16] =   1;
	RoomArray[17] =  -1;
	
	RoomArray[18] =  1;
	RoomArray[19] = -1;
	RoomArray[20] = -1;
	
	RoomArray[21] =  -1;
	RoomArray[22] = -1;
	RoomArray[23] = -1;


	// make norm
	unsigned int col = 0;
	for ( size_t i = 0 ; i < RoomArray.size() ; i=i+3 ){
		VecNSlice<Real> A ( RoomArray.data()+i , 3 );
		Room.push_back( A );
		Real Norm = sqrt( DotProduct( Room[col] , Room[col] ) );
		Room[col] = Room[col] / Norm;
		col++;
	}
}






void AnalyzePotential::AssignVector( std::vector<Real>& data , std::vector<VecNSlice<Real> >& Pointers , unsigned int dim ){
	for ( size_t i = 0 ; i < data.size() ; i = i + dim ){
		VecNSlice<Real> A( data.data() + i , dim );
                Pointers.push_back( A );
	}
}


void AnalyzePotential::PotentialAnaMain( const std::vector<VecNSlice<Real> >& Pos , const std::vector<Real>& Epot , const std::vector<std::vector<Real> >& lattice ,
		                           const std::vector<VecNSlice<Real> >& Force ){


	Delta[0] = get_nearest_image( Center[0] , Pos[0] ) - Center[0];
	for ( size_t i = 1 ; i < Pos.size() ; i++ ){
		std::vector<Real> TempPos =  get_nearest_image( EquiPos[i] , Pos[i] - Delta[0] );
                std::vector<Real> delta   =  TempPos - EquiPos[i];
		// displacement of atom from equilibrium position with removed center of mass drift
                delta = MatVec( lattice , delta );
		Real NormDelta = VecNorm( delta );
                unsigned int bin    =  computeBin( NormDelta , dx , Nbins );
		unsigned int dir = CheckWhichDirection( delta );

		Counter[ dir ][ bin ] ++;
                std::vector<Real> CartForce = MatVec( lattice , Force[i] );
                Real NormForce = VecNorm( CartForce );
                Histograms[dir][bin] +=  Epot[i];
                ForceHistogram[dir][bin] +=  NormForce;
        }
}

unsigned int AnalyzePotential::CheckWhichDirection( const std::vector<Real>& delta ){

	Real Norm = VecNorm( delta );

	std::vector<Real> UnitVec = delta / Norm;

	Real value = -1000;
	unsigned int Dir = 10; 
	for ( size_t i = 0 ; i < Axis.size() ;i++ ){
		Real DotP = DotProduct( Axis[i] , UnitVec );
		if ( DotP > value ){
			value = DotP;
			Dir = 0;
		}
	}
	
	for ( size_t i = 0 ; i < Face.size() ;i++ ){
		Real DotP = DotProduct( Face[i] , UnitVec );
		if ( DotP > value ){
			value = DotP;
			Dir = 1;
		}
	}
	
	for ( size_t i = 0 ; i < Room.size() ;i++ ){
		Real DotP = DotProduct( Room[i] , UnitVec );
		if ( DotP > value ){
			value = DotP;
			Dir = 2;
		}
	}
	return Dir;

}



void AnalyzePotential::PotentialAnaFinalize( void ){


	std::vector<std::vector<Real> > DirDist;
        DirDist.resize( Counter.size() );


        for ( size_t i = 0 ; i < Histograms.size() ; i++ ){
                unsigned int Sum = std::accumulate( Counter[i].begin() , Counter[i].end() , 0 );
                Real norm = D1_simps( std::vector<Real>( Counter[i].begin() , Counter[i].end() ) )*dx;
                DirDist[i].resize( Counter[i].size() );
                for ( size_t j = 0 ; j < Histograms[i].size() ; j++ ){
                        //if ( Counter[i][j] != 0 ) Histograms[i][j] /= Real( Counter[i][j] );
                        Histograms[i][j]  /=  Real( Sum ) / 1000;
                        ForceHistogram[i][j]  /=  Real( Sum );
                        DirDist[i][j]      =  Real( Counter[i][j] ) / norm;

                }
        }

        FILE * outfile;  // energy file
        FILE * outfile2;  // distribution file
        FILE * outfile3;  // force distribution file
        outfile = fopen( "PotentialAnalysis.dat" , "w" );
        outfile2 = fopen( "Distribution.dat" , "w" );
        outfile3 = fopen( "ForceDistribution.dat" , "w" );

        for ( auto i = 0 ; i < Nbins ; i++ ){
             fprintf( outfile , "    %f" , Real( i )*dx+dx*0.5 );
             fprintf( outfile2 , "    %f" , Real( i )*dx+dx*0.5 );
             fprintf( outfile3 , "    %f" , Real( i )*dx+dx*0.5 );
             for ( auto j = 0 ; j < Histograms.size() ; j++ ){
                     fprintf( outfile , "    %f" , Histograms[j][i] );
                     fprintf( outfile2 , "    %f" , DirDist[j][i] );
                     fprintf( outfile3 , "    %f" , ForceHistogram[j][i] );
             }
             fprintf( outfile , "\n" );
             fprintf( outfile2 , "\n" );
             fprintf( outfile3 , "\n" );
        }
        fclose( outfile );
        fclose( outfile2 );
        fclose( outfile3 );
}
