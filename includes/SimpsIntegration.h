#include <vector>


#ifdef USE_DOUBLES
  typedef double Real;
#else
  typedef float Real;
#endif


#ifndef _INTEGRATE
#define _INTEGRATE

 Real D1_simps( std::vector<Real> data );


 Real D2_simps( std::vector<std::vector<Real>> data );
#endif
