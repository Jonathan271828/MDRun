#include <math.h>
#include <vector>
#include <algorithm>
#include <stdio.h>
#include <iostream>


#ifdef USE_DOUBLES
  typedef double Real;
#else
  typedef float Real;
#endif


#ifndef _VectorTools
#define _VectorTools

  template <typename T>
  std::vector<Real> CrossProduct3D( const std::vector<T> A , const std::vector<T> B ){

	  std::vector<T> result;
	  result.resize( 3 );
	  result[ 0 ]  =  A[ 1 ] * B[ 2 ] - A[ 2 ] * B[ 1 ] ;
	  result[ 1 ]  =  A[ 2 ] * B[ 0 ] - A[ 0 ] * B[ 2 ] ;
	  result[ 2 ]  =  A[ 0 ] * B[ 1 ] - A[ 1 ] * B[ 0 ] ;

	  return result;
  }

  template<typename T>
  T Product( const std::vector<T> Vector ){
	  T res = T( 1 );
	  for ( auto i = 0 ; i < Vector.size() ; i++ ){
		  res *= Vector[i];
	  }
	  return res;
  }

  template<typename T>
  std::vector<Real> operator*( const std::vector<T>& x , const T a ){

	  std::vector<T> y ;
	  y.resize( x.size() );
	  for ( auto i = 0 ; i < x.size() ; i++ ){
		  y[i] = x[i] * a;
          }
	  return y;
  }

  template<typename T>
  std::vector<Real> operator*( const T a , const std::vector<T>& x ){

	  return x * a;
  }


  template<typename T>
  std::vector<Real> operator/( const std::vector<T>& x , const T a ){

	  std::vector<T> y ;
	  y.resize( x.size() );
	  for ( auto i = 0 ; i < x.size() ; i++ ){
		  y[i] = x[i] / a;
          }
	  return y;
  }



  template<typename T>
  std::vector<Real> operator+( const std::vector<T>& x , const T a ){

	  std::vector<T> y ;
	  y.resize( x.size() );
	  for ( auto i = 0 ; i < x.size() ; i++ ){
		  y[i] = x[i] + a;
          }
	  return y;
  }

  template<typename T>
  std::vector<Real> operator+( const T a , const std::vector<T>& x ){

	  return x + a;
  }




  // vector + vector
  template<typename T>
  std::vector<Real> operator+( const std::vector<T>& x , const std::vector<T>& a ){

	  std::vector<T> y ;
	  y.resize( x.size() );
	  for ( auto i = 0 ; i < x.size() ; i++ ){
		  y[i] = x[i] + a[i] ;
          }
	  return y;
  }


  template<typename T>
  std::vector<Real> operator-( const std::vector<T>& x , const T a ){

	  std::vector<T> y ;
	  y.resize( x.size() );
	  for ( auto i = 0 ; i < x.size() ; i++ ){
		  y[i] = x[i] - a;
          }
	  return y;
  }

  template<typename T>
  std::vector<Real> operator-( const T a , const std::vector<T>& x ){

	  return x - a;
  }


  // vector - vector
  template<typename T>
  std::vector<Real> operator-( const std::vector<T>& x , const std::vector<T>& a ){

	  std::vector<T> y ;
	  y.resize( x.size() );
	  for ( auto i = 0 ; i < x.size() ; i++ ){
		  y[i] = x[i] - a[i] ;
          }
	  return y;
  }




  // matrix + scalar
  template<typename T>
  std::vector<std::vector<Real>> operator+( const std::vector<std::vector<T>>& x , const T a ){

	  std::vector<std::vector<T>> y ;
	  y.resize( x.size() );
	  for ( auto i = 0 ; i < x.size() ; i++ ){
		  y[i]  = x[i] + a;
          }
	  return y;
  }

  // scalar + matrix
  template<typename T>
  std::vector<std::vector<Real>> operator+( const T a , const std::vector<std::vector<T>> & x ){

	  return x + a;
  }


  // matrix + matrix
  template<typename T>
  std::vector<std::vector<Real>> operator+( const std::vector<std::vector<T>>& x ,
		               const std::vector<std::vector<T>>& a ){

	  std::vector<std::vector<T>> y ;
	  y.resize( x.size() );
	  for ( auto i = 0 ; i < x.size() ; i++ ){
		  y[i]  = x[i] + a[i];
          }
	  return y;
  }


  // matrix - scalar
  template<typename T>
  std::vector<std::vector<Real>> operator-( const std::vector<std::vector<T>>& x , const T a ){

	  std::vector<std::vector<T>> y ;
	  y.resize( x.size() );
	  for ( auto i = 0 ; i < x.size() ; i++ ){
		  y[i]  = x[i] - a;
          }
	  return y;
  }

  // scalar - matrix
  template<typename T>
  std::vector<std::vector<Real>> operator-( const T a , const std::vector<std::vector<T>> & x ){

	  return x - a;
  }


  // matrix - matrix
  template<typename T>
  std::vector<std::vector<Real>> operator-( const std::vector<std::vector<T>>& x ,
		               const std::vector<std::vector<T>>& a ){

	  std::vector<std::vector<T>> y ;
	  y.resize( x.size() );
	  for ( auto i = 0 ; i < x.size() ; i++ ){
		  y[i]  = x[i] + a[i];
          }
	  return y;
  }




  // matrix * scalar
  template<typename T>
  std::vector<std::vector<Real>> operator*( const std::vector<std::vector<T>>& x , const T a ){

	  std::vector<std::vector<T>> y ;
	  y.resize( x.size() );
	  for ( auto i = 0 ; i < x.size() ; i++ ){
		  y[i]  = x[i] * a;
          }
	  return y;
  }

  // scalar * matrix
  template<typename T>
  std::vector<std::vector<Real>> operator*( const T a , const std::vector<std::vector<T>> & x ){

	  return x * a;
  }


  // vector * vector elementwise
  template<typename T>
  std::vector<Real> operator*( const std::vector<T>& x , const std::vector<T>& a ){

	  std::vector<T> y ;
	  y.resize( x.size() );
	  for ( auto i = 0 ; i < x.size() ; i++ ){
		  y[i] = x[i] * a[i] ;
          }
	  return y;
  }


  // matrix * matrix elemtwise multiplication
  template<typename T>
  std::vector<std::vector<Real>> operator*( const std::vector<std::vector<T>>& x ,
		               const std::vector<std::vector<T>>& a ){

	  std::vector<std::vector<T>> y ;
	  y.resize( x.size() );
	  for ( auto i = 0 ; i < x.size() ; i++ ){
		  y[i]  = x[i] * a[i];
          }
	  return y;
  }



 template<typename T>
 std::vector<T> VecAbs( std::vector<T> x ){

	 std::vector<T> y ;
	 y.resize( x.size() );
	 for ( auto i = 0 ; i < x.size() ; i++ ){
		 y[i] = abs( x[i] );
         }
	 return y;
 }




 std::vector<Real> AbsVec( std::vector<Real> x );
  /*
   *
   * Matrix times vector x = A * v
   *
   */

  std::vector<Real> MatVec( std::vector< std::vector<Real >> A , std::vector<Real> v );
  /*
   *
   * scalar product
   *
   */

  Real DotProduct( std::vector<Real> A , std::vector<Real> B );
  /*
   *
   * norm of a vector
   *
   */
  Real VecNorm( std::vector<Real> A );


  /*
   *
   * vector + vector
   *
   */

  std::vector<Real> AddVecs( std::vector<Real> a , std::vector<Real> b );

/*
   *
   * vector - vector
   *
   */

  std::vector<Real> SubVecs( std::vector<Real> a , std::vector<Real> b );

  /*
   *
   * vector times scalar
   *
   */

  std::vector< Real > ScaTVec( Real x , std::vector<Real> a );


  /*
   *
   * vector divide scalar
   *
   */

  std::vector< Real > ScaDivVec( Real x , std::vector<Real> a );


 std::vector<std::vector<Real>> ScaTMat( Real x , std::vector<std::vector<Real>> M );


 //
 // multiply two matrices
 //
 std::vector<std::vector<Real>> MatMul( std::vector<std::vector<Real>> A , std::vector<std::vector<Real>> B );


 //
 // project vector on the x,y or z direction
 //
 std::vector<Real> ProjectXYZ( std::vector<Real> x );
 
 // compute vector exponential

 std::vector<Real> VecExp( std::vector<Real> x );

#endif
