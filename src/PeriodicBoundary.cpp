
#ifndef _PeriodicBound
#include "PeriodicBoundary.h"
#endif 




std::vector<Real> get_nearest_image( std::vector<Real> A , std::vector<Real> B ){

	 if ( A.size() != B.size() ){
	    printf( "Error in function get_nearest_image \n" );
            printf( "sizes of input vectors do not match\n" );
	    return std::vector<Real> ( A.size() , Real(0) );
     }

	 std::vector<Real> C( B.begin() , B.end() );

	 for ( auto i = 0 ; i < A.size() ; i++ ){
		 Real delta = A[i] - B[i];
		 if ( std::abs( delta ) > 0.5 ){
			 C[i] = C[ i ]  +  Sign( delta ) * Real( 1.0 );
		 }
	 }
	 return C;
 }



std::vector<Real> get_nearest_image( VecNSlice<Real> A , VecNSlice<Real> B ){

	 std::vector<Real> C( A.size() );

	 for ( auto i = 0 ; i < A.size() ; i++ ){
		 C[i] = B[i];
		 Real delta = A[i] - B[i];
		 if ( std::abs( delta ) > 0.5 ){
			 C[i] = C[ i ]  +  Sign( delta ) * Real( 1.0 );
		 }
	 }
	 return C;
 }


std::vector<Real> get_nearest_image( std::vector<Real> A , VecNSlice<Real> B ){

	 std::vector<Real> C( A.size() );

	 for ( auto i = 0 ; i < A.size() ; i++ ){
		 C[i] = B[i];
		 Real delta = A[i] - B[i];
		 if ( std::abs( delta ) > 0.5 ){
			 C[i] = C[ i ]  +  Sign( delta ) * Real( 1.0 );
		 }
	 }
	 return C;
 }

std::vector<Real> get_nearest_image( VecNSlice<Real> A , std::vector<Real> B ){

	 std::vector<Real> C( A.size() );

	 for ( auto i = 0 ; i < A.size() ; i++ ){
		 C[i] = B[i];
		 Real delta = A[i] - B[i];
		 if ( std::abs( delta ) > 0.5 ){
			 C[i] = C[ i ]  +  Sign( delta ) * Real( 1.0 );
		 }
	 }
	 return C;
 }



//
// return nearest image of vector B with respect to vector A
// in cartesian coordinates
//
std::vector<Real> get_nearest_image_cart( std::vector<Real>& A , std::vector<Real>& B ,
		                                  const std::vector<std::vector<Real>>& lattice ){

	 if ( A.size() != B.size() ){
	    printf( "Error in function get_nearest_image \n" );
            printf( "sizes of input vectors do not match\n" );
	    return std::vector<Real> ( A.size() , Real(0) );
         }

	 std::vector<Real> tempA = MatVec( lattice , A );
	 Real norm = VecNorm( SubVecs( tempA , MatVec( lattice , B ) ) );

	 std::vector<Real> C = MatVec( lattice , B );

	 std::vector<Real> tempB( B.begin() , B.end() );


	 for ( auto i = -1 ; i < 2 ; i++ ){
		 for ( auto j = -1 ; j < 2 ; j++ ){
			 for ( auto k = -1 ; k < 2 ; k++ ){
				 if ( i==0 && j==0 && k==0 ){continue;}
			         tempB[0] = B[0] + Real( i );
			 	 tempB[1] = B[1] + Real( j );
				 tempB[2] = B[2] + Real( k );
				 Real TempNorm =
					 VecNorm( MatVec( lattice , SubVecs( A , tempB ) ) );
				 if ( TempNorm < norm ){
					 C = tempB;
					 norm = TempNorm;
				 }
			 }
		 }
	 }
	 return C;
 }


 void ApplyPeriodicBound( std::vector<Real>& Pos ){

	 for ( auto i = 0 ; i < Pos.size() ; i++ ){
		 if ( Pos[i] > Real( 1 ) ){
			 Pos[i] = Pos[i] - Real( 1 );
		 }
		 else if( Pos[i] < Real( 0 ) ){
			 Pos[i] = Pos[i] + Real( 1 );
		 }
	 }
}


void ApplyPeriodicBound( VecNSlice<Real>& Pos ){

	 for ( auto i = 0 ; i < Pos.size() ; i++ ){
		 Pos[i] = Pos[i];
		 if ( Pos[i] > Real( 1 ) ){
			 Pos[i] = Pos[i] - Real( 1 );
		 }
		 else if( Pos[i] < Real( 0 ) ){
			 Pos[i] = Pos[i] + Real( 1 );
		 }
	 }
 }


