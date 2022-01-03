#include "MachineLearningFF.h"


MachineLearningForces::MachineLearningForces( std::string file1 , std::string file2 , std::string file3 , 
				                              std::function <void ( std::vector<Real>& data , std::vector<VecNSlice<Real> >& Pointers ,
				   	                          unsigned int dim )> y , 
											  unsigned int Natoms , Real Cutoff ,  
							                  std::vector<VecNSlice<Real> > PosDir , std::vector<VecNSlice<Real> > PosCar ,
							                  std::vector<std::vector<Real> > Box ){


	std::cout << "****************************************" << std::endl;
	std::cout << "Reading machine learning parameters from files" << std::endl; 

    std::cout << "Reading weights from file " << file1 << std::endl;


	// initialize weights
	ReadTable<Real> Reader( file1 );
	Reader.ReadData();
	std::vector<std::vector<Real> > TempData;
	TempData = Reader.ReturnData();
	Weights.resize( TempData.size() );
	for ( auto i = 0 ; i < TempData.size() ; i++ ){
		Weights[ i ] = TempData[i][0] ;
	}


    std::cout << "Reading descriptors from file " << file2 << std::endl;
	// initialize descriptors
	( &Reader )-> ~ReadTable<Real>();   //call the destructor explicitly
	new ( &Reader ) ReadTable<Real>( file2 );
	Reader.ReadData();
	TempData = Reader.ReturnData();
	DescriptorsData.resize( TempData.size() * TempData[0].size() );
    y( DescriptorsData , Descriptors , TempData[0].size() );
	for ( auto i = 0 ; i < TempData.size() ; i++ ){
		for ( auto j = 0 ; j < TempData[0].size() ; j++ ){
			Descriptors[ i ][ j ] = TempData[ i ][ j ];
		}
	}



	NBasis = Descriptors[ 0 ].size();

	Frequencies.resize( NBasis );
	for ( auto i = 0 ; i < NBasis ; i++ ){
		Frequencies[ i ] = M_PI * ( Real( i ) + 1 ) / Cutoff;
	}



    std::cout << "Reading descriptor derivatives from file " << file3 << std::endl;
	// initialize descriptors derivatives
	( &Reader )-> ~ReadTable<Real>();   //call the destructor explicitly
	new ( &Reader ) ReadTable<Real>( file3 );
	Reader.ReadData();
	TempData = Reader.ReturnData();
	DescriptorsDerivativeData.resize( TempData.size() * TempData[0].size() );
    y( DescriptorsDerivativeData , DescriptorsDerivative , TempData[0].size() );

	// make data arrays for descriptors of actual structure
	DescriptorActStrucData.resize( NBasis * Natoms );
    y( DescriptorActStrucData , DescriptorActStruc , NBasis );

	// make data arrays for descriptors derivatives of actual structure
	DescriptorDerivativeActStrucData.resize( NBasis * Natoms * 3 );
    y( DescriptorDerivativeActStrucData , DescriptorDerivativeActStruc , NBasis );


	for ( auto i = 0 ; i < NBasis ; i++ ){
	   FilterFunction.push_back( Cutoff / ( Real( NBasis ) * 1.1 ) * 
	    ( cos( 0.5 * Frequencies[ i ] * 2 * M_PI/ ( Frequencies[ NBasis -1 ] ) ) + 1 ) );
	}




	std::cout << "Machine learning force field initialized" << std::endl;	
	std::cout << "****************************************" << std::endl;

	LinkList.ComputeNNListPeriodicBoundInit( PosCar , PosDir , Cutoff , Box );
	LinkList.PrintNNList();
}




void MachineLearningForces::ComputeDescriptors( std::vector<VecNSlice<Real> >& PosDir , const std::vector<std::vector<Real> >& Box ){


	for ( auto i = 0 ; i < Frequencies.size() ; i++ ){
		for ( auto j = 0 ; j < LinkList.NNList.size() ; j++ ){
			for ( auto k = 0 ; k < LinkList.NNList[ j ].size(); k++ ){
				std::vector<Real> Temp =  get_nearest_image( PosDir[ j ] , PosDir[ LinkList.NNList[ j ][ k ] ] ) - PosDir[ j ] ;
				Temp = MatVec( Box , Temp );
				Real norm = VecNorm( Temp );
				DescriptorActStruc[ j ][ i ] = DescriptorActStruc[ j ][ i ] + sin( Frequencies[i] * norm );
				// compute force descriptors
				for ( auto xyz = 0 ; xyz < 3 ; xyz++ ){
					DescriptorDerivativeActStruc[ j*3 + xyz ][ i ] =  DescriptorDerivativeActStruc[ j*3 + xyz ][ i ] + 
					                                          cos( Frequencies[ i ] * VecNorm( Temp ) ) * 
														      Temp[ xyz ] / norm ;
				}
			}
			DescriptorActStruc[ j ][ i ] = DescriptorActStruc[ j ][ i ] * FilterFunction[ i ];
			for ( auto xyz = 0 ; xyz < 3 ; xyz++ ){
				DescriptorDerivativeActStruc[ j*3 + xyz ][ i ] =
				        	DescriptorDerivativeActStruc[ j*3 + xyz ][ i ] * FilterFunction[ i ];
		    }
		}
	}	
}




void MachineLearningForces::ComputeForces( std::vector<VecNSlice<Real> >& ForceDir , std::vector<VecNSlice<Real> >& ForceCar ,  
		                                   const std::vector<std::vector<Real> > Box ){

	for ( auto j = 0 ; j < ForceDir.size() ; j++ ){
		for ( auto k = 0 ; k < ForceDir[ 0 ].size() ; k++ ){
			ForceCar[ j ][ k ] = Real( 0 );
			// sum over basis set size
			for ( auto i = 0 ; i < Descriptors.size() ; i++ ){
				Real cicj = Real( 0 );
				// Cmpute overlap integral
				for ( auto NB = 0 ; NB < NBasis ; NB++ ){
					cicj =  cicj + Descriptors[ i ][ NB ] * 
						    DescriptorDerivativeActStruc[ j*3 + k ][ NB ];
			    }
		        ForceCar[ j ][ k ] = ForceCar[ j ][ k ] - cicj * Weights[ i ];
		    }
	    }
		ForceDir[ j ] = MatVec( Box , ForceCar[ j ] );
	}
}



Real MachineLearningForces::ComputeEnergy( void ){

	Real energy = Real( 0 );

	// sum over atoms in current structure
	for ( auto j = 0 ; j < DescriptorActStruc.size() ; j++ ){
		// sum over basis set size
		for ( auto i = 0 ; i < Descriptors.size() ; i++ ){
			Real cicj = Real( 0 );
			// Cmpute overlap integral
			for ( auto NB = 0 ; NB < NBasis ; NB++ ){
				cicj =  cicj + Descriptors[ i ][ NB ] * DescriptorActStruc[ j ][ NB ];
			}
		    energy = energy + cicj * Weights[ i ];
		}
	}
	return energy;
}
