

#ifndef _INTEGRATE
#include "SimpsIntegration.h"
#endif



 Real D1_simps( std::vector<Real> data ){

         Real result = 0;
         Real temp1  = 0 ;
         Real temp2 = 0 ;

         for ( auto i = 2 ; i < data.size() ; i+=2 ){
                 temp1 += data[ i ];
                 temp2 += data[ i - 1 ];
         }
         temp1 *= 4.0;
         temp2 *= 2.0;
         result = temp1 + temp2 + data[ data.size() - 1 ] + data[ 0 ];
         result /= 3.0;
         return result;
 }


 Real D2_simps( std::vector<std::vector<Real>> data ){

         std::vector<Real> integral;
	 integral.resize( data.size() );
         Real result;

         for ( auto i = 0 ; i < data.size() ; i++ ){
                 integral[ i ]  =  D1_simps( data[ i ] ) ;
         }
         result = D1_simps( integral );

         return result;
 }

