#include <fstream>
#include <string>
#include <iostream>
#include <string>
#include <sstream>
#include <numeric>
#include <vector>
#include <ios>

#ifdef USE_DOUBLES
  typedef double Real;
#else
    typedef float Real;
#endif





#ifdef USE_DOUBLES
typedef double Real;
#else
typedef float Real;
#endif


#ifndef _FLAGFinder
#define _FLAGFinder


const std::string WHITESPACE = " \n\r\t\f\v";

std::string ltrim(const std::string& s);

std::string rtrim(const std::string& s);

std::string trim(const std::string& s);

Real extractValue( std::string data );

std::vector<Real> extractData( std::string data );
std::string extractStringData( std::string data );



class FlagFinder{


	private:
		std::string fname ;
	    std::ifstream infile;
	public:
		FlagFinder( std::string file );
		~FlagFinder( void ){};
		Real CheckFlag( std::string flag );

		std::vector<Real> CheckFlagVec( std::string flag );

		std::string CheckStringFlag( std::string flag );
 };








#endif 
