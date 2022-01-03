#pragma once

#include <iostream>
#include <vector>
#include <chrono>
#include <math.h>


#include "rand.h"
#ifndef _VectorTools
#include "MatMulVec.h"
#endif
#ifndef _PeriodicBound
#include "PeriodicBoundary.h"
#endif
#ifndef _Forces
#include "force.h"
#endif

#ifndef _MDSPEC
#include "atoms.h"
#endif






#ifdef USE_DOUBLES
typedef double Real;
#else
typedef float Real;
#endif


//
//
//####################################################################
//
//
// this object is completelty fucked up and needs major revision
// don't use it or maybe even delete it
//
//
//####################################################################
//
//


#ifndef _PARTICLE_INSERTION
#define _PARTICLE_INSERTION
class InsertParticle{


	private:
		Ran Rand;
		Real temperature;
		Units unit;
		std::vector<unsigned int> LinkList;
                void DefineLinkList( std::vector<Real> PosNew , std::vector<std::vector<Real> > Pos , Real Cutoff , 
		  		     std::vector<std::vector<Real> > Box );
		unsigned int Nsteps;
		std::vector<Real> NewPosition( unsigned int N );


	public:
		InsertParticle( unsigned int Steps , Real T ):
			Rand( std::chrono::system_clock::to_time_t( std::chrono::system_clock::now() ) )
	                {  
			      Nsteps      =  Steps;
			      temperature =  T * unit.kB; 
			};
		~InsertParticle( void ){};
		std::vector<Real> MonteCarloInsertion( std::vector<VecNSlice<Real> > & Pos ,
			                               std::vector<VecNSlice<Real> > & Velos ,
						       const Real Cutoff , 
				                       const std::vector<std::vector<Real> > Box , 
						       ComputeForces force ,
						       Real Sigma );





};
#endif
