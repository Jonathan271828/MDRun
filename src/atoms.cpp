#ifndef _MDSPEC
#include "atoms.h"
#endif



void AtomType::Transform( std::vector<VecNSlice<Real> >& A , const std::vector<VecNSlice<Real> >& B ,
		                  const std::vector<std::vector<Real> > lattice ){

#ifdef USE_OPENMP
#pragma omp parallel for
#endif
	 for ( auto i = 0 ; i < A.size() ; i++ ){
		 A[i] = MatVec( lattice , B[i] );
     }
 }


 void AtomType::ReadRestartStructure( std::string fname ){


	 std::cout << fname << std::endl;


	 std::ifstream infile( fname.c_str() );

	 std::string data;

	 unsigned int Nline=0;

	 if ( !infile ){
		 std::cout << "############# ERROR #############" << std::endl;
		 std::cout << "      File does not exist        " << std::endl;
		 std::cout << "############# ERROR #############" << std::endl;
		 return;
	 }

	 unsigned int Natom = 0;

	 std::string Temp;
	 while ( std::getline( infile , data ) ){
		 data = trim( data );
		 if ( Nline == 3 ){
		    unsigned int Natoms = ::atof( trim( data ).c_str() );

		    PosDirData.resize( 0 );
		    PosCarData.resize( Natoms * 3 );
		    VeloDirData.resize( Natoms * 3 );
		    VeloCarData.resize( Natoms * 3 );
		    ForceDirData.resize( Natoms * 3 );
		    ForceCarData.resize( Natoms * 3 );
		 }
		 if ( Nline > 8 ){
		    unsigned int col = 0;
                    std::stringstream StrStream( trim( data ).c_str() );
		    while ( StrStream >> Temp ){
		       if ( col > 1 ) PosDirData.push_back( ::atof( Temp.c_str() ) );
		       col ++;
		    }
		    Natom++;
		 }
		 Nline++;
	 }
 }

 void AtomType::ReadStartStruc( std::string fname , SystemParams Params ){

	 std::cout << "Reading start structure from file " << fname << std::endl;


	 std::ifstream infile( fname.c_str() );
	 if ( !infile ){
		 std::cout << "############# ERROR #############" << std::endl;
		 std::cout << "      File does not exist        " << std::endl;
		 std::cout << "############# ERROR #############" << std::endl;
		 return;
	 }




	 std::string data;
	 std::string Temp;
	 unsigned int Nline = 0;
	 while ( std::getline( infile , data ) ){
		 data = trim( data );
		 std::stringstream StrStream( trim( data ).c_str() );

		 if ( Nline > 1 && Nline < 5 ){
			 unsigned int counter = 0;
			 while ( StrStream >> Temp ){
                Params.Box[ Nline - 2 ][ counter ] = ::atof( Temp.c_str() );
                Params.RecBox[ Nline - 2 ][ counter ]  =  Real( 1 ) / Params.Box[ Nline - 2 ][ counter ] ;
				counter ++;
			 }
		 }



		 if ( Nline == 6 ){
		    unsigned int Natoms = ::atof( trim( data ).c_str() );
		    PosDirData.resize( 0 );
		 }
		 if ( Nline > 7 ){
                    std::stringstream StrStream( trim( data ).c_str() );
		    while ( StrStream >> Temp ){
		       PosDirData.push_back( ::atof( Temp.c_str() ) );
		    }
		 }
		 Nline++;
	 }



	 for ( auto i = 0 ; i < PosDir.size() ; i++ ){
		 PosCar[ i ] = MatVec( Params.Box , PosDir[ i ] );
	 }

	 for ( auto i = 0 ; i < Params.Box.size() ; i++ ){
		 for ( auto j = 0 ; j < Params.Box.size() ; j++ ){
			 if ( i != j ){
				 Params.RecBox[ i ][ j ] = Real( 0 );
			 }
		 }
	 }
 }



 void AtomType::InitializePointers( std::vector<Real>& data , std::vector<VecNSlice<Real> >& Pointers ,
		                    unsigned int dim ){
	 //Pointers.resize( 0 );
	 for ( auto i = 0 ; i < data.size() ; i=i+dim ){
		 VecNSlice<Real> A( data.data() + i , dim );
		 Pointers.push_back( A );
	 }
 }


HeatGradientStruc & HeatGradientStruc::operator=( const HeatGradientStruc & InData ){    
                                                        
      if ( this != & InData ){
	      this -> HeatGradient   =  InData.HeatGradient;
              this -> temperatureA   =  InData.temperatureA;
              this -> temperatureB   =  InData.temperatureB;
              this -> Direction      =  InData.Direction;
	      this -> nu             =  InData.nu;
	      this -> ActiveMeasure  =  InData.ActiveMeasure;
	      this -> EquiSteps      =  InData.EquiSteps;

      }
      return *this;
}                                                            
