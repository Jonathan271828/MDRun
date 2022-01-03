#ifndef _VectorTools
#include "MatMulVec.h"
#endif

 std::vector<Real> AbsVec( std::vector<Real> x ){

	 std::vector<Real> y ;
	 y.resize( x.size() );
	 for ( auto i = 0 ; i < x.size() ; i++ ){
		 y[i] = abs( x[i] );
	 }
	 return y;
 }

  /*
   *
   * Matrix times vector x = A * v
   *
   */

  std::vector<Real> MatVec( std::vector< std::vector<Real > > A , std::vector<Real> v ){

	  std::vector<Real> Result;
	  Result.resize( v.size() );
	  for ( auto i = 0 ; i < A.size() ; i++ ){
		  for ( auto j = 0 ; j < A[i].size() ; j++ ){
			  Result[ i ]  =  Result[ i ] + A[i][j] * v[j];
		  }
	  }

	  return Result;
  }


  /*
   *
   * scalar product
   *
   */

  Real DotProduct( std::vector<Real> A , std::vector<Real> B ){

	  Real dotpi = Real( 0 );
	  for ( auto i = 0 ; i < A.size() ; i++ ){
		  dotpi = dotpi + A[i] * B[i] ;
	  }
	  return dotpi;
  }

  /*
   *
   * norm of a vector
   *
   */
  Real VecNorm( std::vector<Real> A ){
	  Real norm = Real( 0 );
	  for ( auto i = 0 ; i < A.size() ; i++ ){
		  norm = norm + A[i] * A[i] ;
	  }
	  return sqrt( norm );
  }


  /*
   *
   * vector + vector
   *
   */

  std::vector<Real> AddVecs( std::vector<Real> a , std::vector<Real> b ){

	  std::vector<Real> c ;

	  for ( auto i = 0 ; i < a.size() ; i++ ){
		  c.push_back( a[ i ] + b[ i ] );
	  }
	  return  c;
  }

/*
   *
   * vector - vector
   *
   */

  std::vector<Real> SubVecs( std::vector<Real> a , std::vector<Real> b ){

	  std::vector<Real> c ;

	  for ( auto i = 0 ; i < a.size() ; i++ ){
		  c.push_back( a[ i ] - b[ i ] );
	  }

	  return  c;
  }

  /*
   *
   * vector times scalar
   *
   */

  std::vector< Real > ScaTVec( Real x , std::vector<Real> a ){

	  std::vector<Real> b;
	  for ( auto i = 0 ; i < a.size() ; i++ ){
		  b.push_back( x * a[ i ] );
	  }
	  return b;
  }


  /*
   *
   * vector divide scalar
   *
   */

  std::vector< Real > ScaDivVec( Real x , std::vector<Real> a ){

	  std::vector<Real> b;
	  for ( auto i = 0 ; i < a.size() ; i++ ){
		  b.push_back( a[ i ] / x );
	  }
	  return b;
  }


 std::vector<std::vector<Real>> ScaTMat( Real x , std::vector<std::vector<Real>> M ){

	 std::vector<std::vector<Real>> A;
	 A.resize( M.size() );
	 for ( auto i = 0 ; i < M.size() ; i++ ){
		 A[i].resize( M[i].size() );
		 for ( auto j = 0 ; j < M[i].size() ; j++ ){
			 A[i][j] = x * M[i][j];
		 }
	 }
	 return A;
 }



 //
 // multiply two matrices
 //
 std::vector<std::vector<Real>> MatMul( std::vector<std::vector<Real>> A , std::vector<std::vector<Real>> B ){

	 if ( A[0].size() != B.size() ){
		 printf( "               ERROR      \n" );
		 printf( "Error in function MatMul dimension of arrays don't fit \n" );
		 printf( "               ERROR      \n" );
		 throw;
	 }

	 std::vector<std::vector<Real>> C;
	 C.resize( A.size() );
	 for ( auto i = 0 ; i < A.size() ; i++ ){
		 C[i].resize( B[ 0 ].size() );
		 for ( auto j = 0 ; j < B.size(); j++ ){
			 C[i][j] = Real( 0 );
			 for ( auto k = 0 ; k < A[i].size() ; k++ ){
				 C[i][j]  =  C[i][j]  +  A[i][k] * B[k][j];
			 }
		 }
	 }
	 return C;
 }

 //
 // project vector on the x,y or z direction
 //
 std::vector<Real> ProjectXYZ( std::vector<Real> x ){

	 std::vector<Real> tmp;
	 std::vector<Real> y;

	 tmp = AbsVec( x ) / VecNorm( x );

	 int Max = std::max_element( tmp.begin() , tmp.end() ) - tmp.begin();
	 if ( Max == 0 ){
		 y = { 1 , 0 , 0 };
	 }
	 else if ( Max == 1 ){
		 y = { 0 , 1 , 0 };
	 }
	 else if ( Max == 2 ){
		 y = { 0 , 0 , 1 };
	 }

	 return y;
 }


 // compute vector exponential

 std::vector<Real> VecExp( std::vector<Real> x ){

	 std::vector<Real> y;
	 y.resize( x.size() );
	 for ( auto i = 0 ; i < x.size() ; i++ ){
		 y[i] = exp( x[i] );
	 }
	 return y;
 }
