#pragma once
#include <math.h>
#include <iostream>

typedef unsigned long long int Ullong;




#ifdef USE_DOUBLES
typedef double Real;
#else
typedef float Real;
#endif




/*****************************************
 *
 * initilize with some random number as seed
 * int64 can be called to return a 
 * 64bit integer random number
 * int32 return a 32 bit random number
 * doub return a double real random
 * number in interval 0,1
 * the period of the generator is roughly
 * ~3.138*10**57
 * programmed after
 * Numerical Recipes for cpp
 * The art of Scientific computing 
 * third edition
 * Cambridge university press
 *
 ****************************************/
 class Ran{

	private:
		Ullong u;
		Ullong v;
		Ullong w;
	public:
		// constructoor call with a Ullong integer as seed
		Ran( Ullong j ): v( 4101942887655102017LL ) , w( 1 ){
			u = j ^ v ;
			int64();
			v = u;
			int64();
			w = v ;
			int64();
		};
		// all the other generators are based on this one
		// creating  64int random number
		inline Ullong int64(){
			u  = u * 2862933555777941757LL + 7046029254386353087LL;
			v ^= v >> 17;
			v ^= v << 31;
			v ^= v >> 8;
			w = 429495766U * ( w & 0xffffffff ) + ( w >> 32 );
			Ullong x = u ^ ( u << 21 );
			x ^= x >> 35;
			x ^= x << 4;
			return ( x + v ) ^ w ;

		};

		// generate uniform random distribution
		inline double doub(){ return 5.42101086242752217E-20 * int64() ;};
		// 32int random number generator
		inline unsigned int int32() { return ( unsigned int ) int64() ;};
		// generate gaussian random distribution
                inline double generate_gauss_dist( double mu , double sigma ){
			double u1 ;
	                double u2 ;
                        double z1 ;
	                static double z2 ;
	                double s ;
	                thread_local bool generate_new ;


	                generate_new  =  !generate_new ;

                        //sigma *= ( double ) 2549.876663134 ;

	                if ( !generate_new ){ 
	                    return z2 * sigma + mu ; }

	                do {
	                    u1  =  doub() * 2.0 - 1.0 ;
	                    u2  =  doub() * 2.0 - 1.0 ;
	                    s = u1 * u1 + u2 * u2 ;
                        } while ( ( s >= 1.0 ) || ( s == 0.0 ) ) ; 
    

	                s = sqrt( -2.0 * log( s ) / s ) ; 

	                z1  =  u1 * s ;
	                z2  =  u2 * s ;

	                return z1  *  sigma +  mu ;
		}
                
		inline double generate_exp_dist(  double lambda ){


			double u;
			u  =  doub();
			return -log( 1.0 - u ) / lambda;
		}


		// draw a gamma distributed random varibale
		// methods is implemented according to
		// George Marsagli and Wai Wan Tsang
		// A Simple Method for Generating Gamma
                // Variables, 
		// ACM Transactions on Mathematical Software,
		// Vol. 26, No. 3, September 2000, Pages 363–372

		inline double generate_gamma_dist( double a ){
			double c;
			double d;
			double x;
			double v;
			double u;
			double value;

			d  =  a - 1.0/3.0;
			c  =  1.0 / sqrt( 9.0 * d );

			// initialize
			u  =   doub();
			x  =   generate_gauss_dist( 0.0 , 1.0 );
                        v  =   computeVGamma( c , x );
			do {
				x  =  generate_gauss_dist( 0.0 , 1.0 );
				v  =  computeVGamma( c , x );
                        } while ( v <= 0 );

			value  =  v * d;
			do {
				u = doub();

				do {
				   x  =  generate_gauss_dist( 0.0 , 1.0 );
				   v  =  computeVGamma( c , x );
				} while ( v <= 0 );

				//std::cout << computeTerm1Gamma( x ) << "   "
				//	  << computeTerm2Gamma( d , v , x ) << "   "
				//	  << log( u  ) << std::endl;
				value = v * d;
			} while ( u > computeTerm1Gamma( x ) ||
				  log( u  ) > computeTerm2Gamma( d , v , x ) );
			return value;
		}

		/// gamma distribution helper routines
		double computeVGamma( double c , double x ){
			double results = 1.0 + c*x;
			return results * results * results; 
		}
		
		double computeTerm1Gamma( double x ){
			return 1.0 - 0.331 * ( x * x )*( x * x );
		}

		double computeTerm2Gamma( double d , double v , double x ){
			return 0.5 * ( x * x ) + d  - d * v + d*log( v );
		}
		// end generating gamma distributed random numbers
};


















