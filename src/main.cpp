#include <iostream>
#include <stdio.h>
#include <string>
#include <numeric>
#include <fstream>
#include <iomanip>

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



#ifndef _MORSESOLID
#include "mainVariants.h"
#endif


int main( int argc , char *argv[] ){


	unsigned int SimulationType; 


	if ( argc < 2 ){
		ErrorMessage( "No input file supplied\n Please supply input file as first command line argument" ,
			          "Main-Initialize" );
		return -1;
	}
	else{
	    FlagFinder Find( argv[1] );
            SimulationType = int( Find.CheckFlag( "Mode" ) );
        }
		
	MainStyles Simulation;

	if ( SimulationType == 1 || SimulationType == 2 ){
	    int status = Simulation.MorseSolidMain( argv[1] , SimulationType );
	}
	if ( SimulationType == 3 ){
	    int status = Simulation.MachineLearningFF( argv[ 1 ] );
	}
	if ( SimulationType == 4 ){
	    int status = Simulation.HarmonicDoubleWellPerturbation( argv[1] );
	}

    return 0;

}
