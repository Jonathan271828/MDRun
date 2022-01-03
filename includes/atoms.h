#include <vector>
#include <math.h>
#include <iostream>
#include <string>
#include <sstream>


#ifndef _VectorTools
#include "MatMulVec.h"
#endif

#ifndef _FLAGFinder
#include "ReadFlags.h"
#endif



#ifdef USE_DOUBLES
typedef double Real;
#else
typedef float Real;
#endif



#ifndef _MDSPEC
#define _MDSPEC




template<typename T>
class VecNSlice
{

    private:
        T * VecN;
        size_t N;

    public:
        VecNSlice( T * address , size_t N ):N( N ) , VecN( address ){};
	VecNSlice( void ){};
    
    T operator []( size_t i ) const {
	    if( i < N ){ 
		    return VecN[ i ];
	    }else{ 
	            std::cout << "Error: Size of vector too large" << std::endl;
                    throw;
	    } ;
    }
    T & operator []( size_t i ) { 
	    if ( i < N ){ 
		    return VecN[ i ]; 
	    }else{ 
		    std::cout << "Error: Size of vector too large" << std::endl;
                    throw; 
	    } ;
    }

    VecNSlice & operator=( const std::vector<T>& input ) 
    {
    	        for ( auto i = 0 ; i < N ; i++ ){
			VecN[i] = input[i];
    		}
		return *this;
    }



    VecNSlice & operator=( const VecNSlice<T>& input ){
            for ( auto i = 0 ; i < N ; i++ ){
                VecN[i] = input[i];
            }
	    return *this;
    }
    

    VecNSlice & operator=( const T input ){
            for ( auto i = 0 ; i < N ; i++ ){
                VecN[i] = input;
            }
	    return *this;
    }
   
   
    std::vector<T> operator+( VecNSlice<T> input )
    {
        std::vector<T> result( N );
            for ( auto i = 0 ; i < N ; i++ ){
                result[i] = VecN[i] + input[i];
            }
        return result;
    }
   
    std::vector<T> operator * ( T scalar ){
        std::vector<T> result( N );
        for ( auto i = 0 ; i < N ; i++ ){
            result[ i ] = VecN[ i ] * scalar;
        }
        return result;
    }
    size_t size( void ) const { return N;}

    void Print( void )const{
	    for ( auto i = 0 ; i < N ; i++ ){
	         std::cout << 
			   this -> VecN[i] << "  ";
	    }
	    std::cout << std::endl;
    }
};

 



//#pragma omp declare reduction( VecNSlicePlus : VecNSlice<Real> : \
//        std::transform( omp_out.begin(), omp_out.end(), omp_in.begin(), omp_out.begin(), std::plus<Real>())) \
//            initializer(omp_priv = decltype(omp_orig)(omp_orig.size()))




template <typename T>
std::vector<T> MatVec( const std::vector<std::vector<T> > & A , const VecNSlice<T> & v ){

	std::vector<T> result( A.size() );
	for ( auto i = 0 ; i < A.size() ; i++ ){
		for ( auto j = 0 ; j < v.size() ; j++ ){
			result[ i ]  =  result[ i ] + A[i][j] * v[j];
		}
	}
	return result;
}


// dot product vecNSlice * VecNSlice
template<typename T>
T DotProduct( VecNSlice<T> A , VecNSlice<T> B ){

	T result = T(0);
	for ( auto i = 0 ; i < B.size() ; i++ ){
		result = result + A[i] * B[i];
	}
	return result;
}

// dot product vecNSlice * std::vector
template<typename T>
T DotProduct( VecNSlice<T> A , std::vector<T> B ){

	T result = T(0);
	for ( auto i = 0 ; i < B.size() ; i++ ){
		result = result + A[i] * B[i];
	}
	return result;
}


// dot product std::vector * VecNSlice
template<typename T>
T DotProduct( std::vector<T> B , VecNSlice<T> A ){

	T result = T(0);
	for ( auto i = 0 ; i < B.size() ; i++ ){
		result = result + A[i] * B[i];
	}
	return result;
}




 // VecNSlice - VecNSlice



 // VecNSlice - VecNSlice
 template<typename T>
 std::vector<T> operator-( VecNSlice<T> x , VecNSlice<T> a ){

     std::vector<T> y ;
     y.resize( x.size() );
     for ( auto i = 0 ; i < x.size() ; i++ ){
         y[i] = x[i] - a[i] ;
     }
     return y;
 }

 // vector - VecNSlice
 template<typename T>
 std::vector<T> operator-( const std::vector<T> x , VecNSlice<T> a ){

     std::vector<T> y ;
	 y.resize( x.size() );

     for ( auto i = 0 ; i < x.size() ; i++ ){
         y[i] = x[i] - a[i] ;
     }
     return y;
 }


 // VecNSlice - vector
 template<typename T>
 std::vector<T> operator-( VecNSlice<T> x , const std::vector<T> a ){

     std::vector<T> y ;
     y.resize( x.size() );
     for ( auto i = 0 ; i < x.size() ; i++ ){
         y[i] = x[i] - a[i] ;
     }
     return y;
 }


 // scalar * VecNSlice
 template<typename T>
 std::vector<T> operator*( T scalar , VecNSlice<T> x ){

     std::vector<T> y ;
     y.resize( x.size() );
     for ( auto i = 0 ; i < x.size() ; i++ ){
         y[i] = x[i] * scalar ;
     }
     return y;
 }

 // VecNSlice * scalar
 template<typename T>
 std::vector<T> operator*( VecNSlice<T> x , T scalar ){

     std::vector<T> y ;
     y.resize( x.size() );
     for ( auto i = 0 ; i < x.size() ; i++ ){
         y[i] = x[i] * scalar ;
     }
     return y;
 }

 // VecNSlice / scalar
 template<typename T>
 std::vector<T> operator/( VecNSlice<T> x , T scalar ){

     std::vector<T> y ;
     y.resize( x.size() );
     for ( auto i = 0 ; i < x.size() ; i++ ){
         y[i] = x[i] / scalar ;
     }
     return y;
 }


 // VecNSlice + vector
 template<typename T>
 std::vector<T> operator+( VecNSlice<T> x , const std::vector<T> a ){

     std::vector<T> y ;
     y.resize( x.size() );
     for ( auto i = 0 ; i < x.size() ; i++ ){
         y[i] = x[i] + a[i] ;
     }
     return y;
 }


 // VecNSlice - vector
 template<typename T>
 std::vector<T> operator+( const std::vector<T> a , VecNSlice<T> x ){

     std::vector<T> y ;
     y.resize( x.size() );
     for ( auto i = 0 ; i < x.size() ; i++ ){
         y[i] = x[i] + a[i] ;
     }
     return y;
 }

 template<typename T>
 T VecNSliceNorm( const VecNSlice<T> & x ){

	 T norm = T(0);
	 for ( auto i = 0 ; i < x.size() ; i++ ){
		 norm  =  norm + x[i] * x[i];
	 }

	 return sqrt( norm );
 }



template<typename T>
void AssociateVecNSlice( std::vector<T>& data ,
                         std::vector<VecNSlice<T> >& Pointers ,
			 unsigned int dim ){

	for ( auto i = 0 ; i < data.size() ; i=i+dim ){
		VecNSlice<Real> A( data.data() + i , dim );
                Pointers.push_back( A );
        }
}





struct HeatGradientStruc{

	      int HeatGradient;      // check if heat gradient is applied
	      Real temperatureA;      // heat gradient temperature A
	      Real temperatureB;      // heat gradient temperature B
	      int Direction; // Direction to which the heat gradient is applied
	      Real nu;
	      int ActiveMeasure; // check if active measurement is switched on
	      int EquiSteps;     // number of steps for equilibration
              HeatGradientStruc & operator=( const HeatGradientStruc & InData );
};





struct SystemParams
      {
	      std::vector<std::vector<Real>> Box;     // Box
	      std::vector<std::vector<Real>> RecBox;  // reciprocal box
	      Real temperature;             // temperature for Canonical Ens
	      Real timeStep;                // time for MD
	      std::vector<int> NParticles;  // number particles x
	      int dimension;                // dimension
	      std::vector<Real> deltaX ;    // Grid spacing
	      int Natoms ;                  // Number of atoms
	      unsigned int Nsteps;          // md steps
	      unsigned int AnaFrequ=100;    // number of

	      Real Tact;      // actual temperature
	      Real Epot;      // actual potential energy
	      Real Ekin;      // actual kinetic energy
	      Real Etot;      // actual total energy
              Real Pressure;  // actual pressure
   	      Real AndersenFrequ ; // collision probability andersen thermo
	      std::string RestartFile;
	      unsigned int Thermostyle; // type of thermostat

	      Real NoseMass ; //Nose-Hoover virtual mass
	      Real LangevinFriction; // friction coefficient Langevin thermostat

	      Real BDPThermoTau;  // chracteristic time of Bussi Donadino Parrinello thermostat

	      HeatGradientStruc HeatVars; // variables needed for heat transport measurement 

              std::vector<Real> DeltaPhonon;
              bool PhononOnOff;         // define if harmonic phonon stuff should be computed
	      int Output;      // defines used output type for trajectories
	      int StartSample; // determines how much MD steps are done before starting sampling
	      unsigned int SwitchThermo;// switches to microcannonical after SwitchThermo steps 
	      std::vector<Real> NNMorseParams; // potential parameters for nearest neighbor morse
	      Real Cutoff;
	      int ParticleInsert ; // insert particle every N steps
	      bool ParticleInsertOnOff ;    // insert particle every N steps
	      
	      std::vector<Real> NNHarmonicParams;     // potential parameters for nearest neighbor morse
	      std::vector<Real> ONHarmonicParams;     // potential parameters for on site harmonic pot
	      std::vector<Real> NNDoubleWellParams;   // potential paramaters for double well potential
	      Real ForceMixxingParameter;             // mixing parameter [0,1] to mix the harmonic and double well potential
};











 struct AtomType
      {


                  std::vector<VecNSlice<Real> > PosDir ;   //position in direct coords
                  std::vector<VecNSlice<Real> > PosCar ;   //position in cartesian coords
		  std::vector<VecNSlice<Real> > VeloDir;   //velocity direct coordinates
		  std::vector<VecNSlice<Real> > VeloCar;   //velocity cartesian cooordinates
		  std::vector<VecNSlice<Real> > ForceDir;   //force direct cooordinates
		  std::vector<VecNSlice<Real> > ForceCar;   //force cartesian cooordinates



                  std::vector<Real> PosDirData ;   //position in direct coords 3N
                  std::vector<Real> PosCarData ;   //position in cartesian coords 3N
		  std::vector<Real> VeloDirData;   //velocity direct coordinates 3N
		  std::vector<Real> VeloCarData;   //velocity cartesian cooordinates 3N
		  std::vector<Real> ForceDirData;   //force direct cooordinates 3N
		  std::vector<Real> ForceCarData;   //force cartesian cooordinates 3N




		  std::vector<Real> Epot ;  // potential energy per particle
		  std::vector<Real> Ekin ;  // kinetic energy per particle
		  std::vector<Real> Mass;         // mass of considered atom
		  std::vector<int> Type;          // integer type indicator
		  std::vector<std::string> Spec;  // atom specification


		  void Transform( std::vector<VecNSlice<Real> >& A , 
				  const std::vector<VecNSlice<Real> >& B ,
    	                          const std::vector<std::vector<Real>> lattice );
		  void ReadRestartStructure( std::string fname );
                  void ReadStartStruc( std::string fname , SystemParams Params );
		  void InitializePointers( std::vector<Real>& data , 
				           std::vector<VecNSlice<Real> >& Pointers ,
   	                                   unsigned int dim );

 };







/*
 *
 * the units used in this molecular dynamics package
 * are angstroem for distance eV for the energy
 * K for the temperature
 * the time step is exotic ( code seconds ;-)= ), but on input and output it
 * will be pico seconds
 * if you want to print the output in ps divide time by timeConv
 * the mass used in this code are Dalton what is eqivalent to the numerical
 * value of molar mass [g/mol] the remaining units are derived from those
 *
 * it is nice to know that these are the same units as used by lammps
 * in the metals style
 */

struct Units
       {
	       Real kB  =  Real( 8.6173303) / Real( 100000 );    //Boltzmann constant
	       // divide picosecond input by this number to
	       // obtain the proper units the code is using
	       // at output just multiply by this factor
	       Real timeConv = 98.22694538739161;   // time conversion; code uses 1cs = 98.22694538739161ps
       };

#endif
