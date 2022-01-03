

#ifndef _FLAGFinder
#include "ReadFlags.h"
#endif



std::string ltrim(const std::string& s)
{
	        size_t start = s.find_first_not_of( WHITESPACE );
			            return ( start == std::string::npos ) ? "" : s.substr(start);
}

std::string rtrim(const std::string& s)
{
	        size_t end = s.find_last_not_of( WHITESPACE );
			            return ( end == std::string::npos ) ? "" : s.substr( 0 , end + 1 );
}

std::string trim(const std::string& s)
{
	        return rtrim(ltrim(s));
}


 Real extractValue( std::string data ){

	 Real value ;

	 size_t found = data.find( "=" ); 
	 if (found != std::string::npos ){
		 value = ::atof( ( trim( data.substr( found + 1 , data.length() ) ) ).c_str() );
	 }
	 return value; 
 }

 std::vector<Real> extractData( std::string data ){
     std::string Temp;
     size_t found = data.find( "=" ); 
     std::vector<Real> RealData ;

     if ( found != std::string::npos ){
	std::stringstream StrStream( trim( data.substr( found + 1 , data.length() ) ).c_str() );
        while ( StrStream >> Temp ){
            RealData.push_back( ::atof( Temp.c_str() ) );
        } 
     }
    return RealData;
 }


 std::string extractStringData( std::string data ){

     std::string Temp;
     size_t found = data.find( "=" ); 
     std::string Result;

     if ( found != std::string::npos ){
	 Result = trim( data.substr( found + 1 , data.length() ) );
     }

     return Result;
 }





 FlagFinder::FlagFinder( std::string file ){

	 fname = file;
	 infile.open ( fname.c_str(), std::ifstream::in);

 }

 
 Real FlagFinder::CheckFlag( std::string flag ){


	std::string data;
	Real x;

	unsigned int len = flag.length();


	infile.clear();
	infile.seekg(0, std::ios::beg);
	while ( std::getline( infile , data ) ){
		data = trim( data );
		if ( data.substr( 0, len ).compare( flag ) == 0 ) {
			x = extractValue( data );
	        infile.clear();
	        infile.seekg(0, std::ios::beg);
			return x;
		}
	}
	std::cout << "########### WARNING ############" << std::endl;
	std::cout << "problem occoured in flag finder" << std::endl; 
	std::cout << "flag " << flag << " not found" << std::endl; 
	std::cout << "########### WARNING ############" << std::endl;
        return Real( -10000000 );
 }
 
 std::vector<Real> FlagFinder::CheckFlagVec( std::string flag ){


	std::string data;
	std::vector<Real> x;

	unsigned int len = flag.length();


	infile.clear();
	infile.seekg(0, std::ios::beg);
	while ( std::getline( infile , data ) ){
		data = trim( data );
		if ( data.substr( 0, len ).compare( flag ) == 0 ) {
		    x = extractData( data );
	            infile.clear();
	            infile.seekg(0, std::ios::beg);
		    return x;
	        }
	}
	std::cout << "########### WARNING ############" << std::endl;
	std::cout << "problem occoured in flag finder" << std::endl; 
	std::cout << "flag " << flag << " not found" << std::endl; 
	std::cout << "########### WARNING ############" << std::endl;
	return x;
 }
 


 std::string FlagFinder::CheckStringFlag( std::string flag ){


	std::string data;
	std::string x;

	unsigned int len = flag.length();


	infile.clear();
	infile.seekg(0, std::ios::beg);
	while ( std::getline( infile , data ) ){
		data = trim( data );
		if ( data.substr( 0, len ).compare( flag ) == 0 ) {
		     x = extractStringData( data );
	             infile.clear();
	             infile.seekg(0, std::ios::beg);
		     return x;
		}
	}
	std::cout << "########### WARNING ############" << std::endl;
	std::cout << "problem occoured in flag finder" << std::endl; 
	std::cout << "flag " << flag << " not found" << std::endl; 
	std::cout << "########### WARNING ############" << std::endl;
    return "NOTFOUND"; 
 }
