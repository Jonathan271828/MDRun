#include <iostream>
#include <stdio.h>
#include <string>
#include <numeric>
#include <fstream>
#include <iomanip>
#include <math.h>
#include <ctime>
#include <chrono>
#include <omp.h>


#include "rand.h"
#include "MachineLearningFF.h"


#ifndef _MDSPEC
#include "atoms.h"
#endif
#ifndef _VectorTools
#include "MatMulVec.h"
#endif
#ifndef _Forces
#include "force.h"
#endif
#ifndef _ANDERSEN
#include "thermostats.h"
#endif
#ifndef _OUTPUT
#include "output.h"
#endif
#ifndef _ANALYZE
#include "analyze.h"
#endif

#ifndef _FLAGFinder
#include "ReadFlags.h"
#endif

#ifndef _Phononsfidi
#include "PhononFD.h"
#endif


#ifndef _PARTICLE_INSERTION 
#include "ParticleInsertion.h"
#endif







#ifdef USE_DOUBLES
typedef double Real;
#else
typedef float Real;
#endif






#ifndef _MORSESOLID
#define _MORSESOLID
class MainStyles
   {
	   public:
		   MainStyles( void ){};
                   int MorseSolidMain( const char * fname , const unsigned int SymType );
		   int MachineLearningFF( const char * fname );
                   int HarmonicDoubleWellPerturbation( const char * fname );
	   private:
		   void SetUpSimulation( AtomType & atoms , SystemParams & System , std::string inputFile );
		   void SetUpHarmonicDoubleWellPerturbation( AtomType & atoms , 
				                        SystemParams & System , 
							std::string inputFile );
		   int MorseListUpdate( AtomType & atoms , SystemParams & System );
};
#endif
