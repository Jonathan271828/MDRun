#include <math.h>
#include <vector>

#include <iostream>
#include <stdio.h>
#include <chrono>
#include <omp.h>
#include <bits/stdc++.h> 


#include "rand.h"
#include "MachineLearningFF.h"


#ifndef _VectorTools
#include "MatMulVec.h"
#endif
#ifndef _MDSPEC
#include "atoms.h"
#endif
#ifndef _Forces
#include "force.h"
#endif
#ifndef _LinkList
#include "LinkList.h"
#endif





#ifdef USE_DOUBLES
typedef double Real;
#else
typedef float Real;
#endif


template<typename T>
std::vector<std::pair<T,int> > sortArr( const std::vector<T> & v ){ 
  
        // Vector to store element 
        // with respective present index 
	std::vector<std::pair<T, int> > vp; 
  
    // Inserting element in pair vector 
    // to keep track of previous indexes 
    for ( auto i = 0; i < v.size() ; ++i ){ 
        vp.push_back( std::make_pair( v[i] , i ) ); 
    } 
  
    // Sorting pair vector 
    sort(vp.begin(), vp.end()); 
    return vp;
} 






 /*
  *
  *
  *
  *
  * perform an Andersen step
  * under the use of a velocity verlet algorithm
  *
  *
  *
  *
  *
  */

#ifndef _ANDERSEN
#define _ANDERSEN
class ThermostatsVeloVerlet {


	public:
		FILE *DebugFile;
		ThermostatsVeloVerlet( void ):
			// constructor random number generator
			Rand( std::chrono::system_clock::to_time_t( std::chrono::system_clock::now() ) )
	                { };


		Units units;
		//ThermostatsVeloVerlet( void ){};
		void InitAndersen( const Real temp , const Real nu , const std::vector<Real>& mass );
		void InitNoseHoover( Real q , unsigned int DOF , Real temp , const std::vector<Real>& mass );
		void InitLangevin( Real Fric , Real temp , Real dt , unsigned int dim ,
				   const std::vector<Real>& Mass );
		void InitHeatGradientThermoAndersen( const HeatGradientStruc& InData , 
					             const std::vector<int>& Ndirs ,
					             const std::vector<Real>& masses );
                void InitHeatGradientThermoLangevin( const HeatGradientStruc& InData ,
		                                     Real Fric , Real dt , unsigned int dim , 
			                             const std::vector<int>& NDirs , 
						     const std::vector<Real>& Mass );
		void InitBDP( Real tau , Real Temp , Real dt , int DOF ,
				const std::vector<Real>& mass );
                void InitHeatGradientThermoBDP( const HeatGradientStruc& InData , Real tau , Real dt , 
		                                const std::vector<Real>& mass , const std::vector<int>& NDirs,
						int dim );



                void InitVelocities( std::vector<VecNSlice<Real> >& Velos , int dim ,
				             const std::vector<std::vector<Real> >& Box ,
				             const std::vector<Real>& Mass );

                void ThermostatIntegrate( std::vector<VecNSlice<Real> >& Pos ,
                                          std::vector<VecNSlice<Real> >& Force ,
                                          std::vector<VecNSlice<Real> >& Velo ,
                                          Real dt ,
                                          const std::vector<std::vector<Real> >& Box , 
                 		          const std::vector<std::vector<Real> >& RecBox ,
					  ComputeForces GetForce , 
		         	          const std::vector<Real>& Mass );
                

		void ThermostatIntegrateHarmonicDoubleWellPerturbation( 
				          std::vector<VecNSlice<Real> >& Pos ,
                                          std::vector<VecNSlice<Real> >& Force ,
                                          std::vector<VecNSlice<Real> >& Velo ,
                                          Real dt ,
                                          const std::vector<std::vector<Real> >& Box , 
                 		          const std::vector<std::vector<Real> >& RecBox ,
					  ComputeForces GetForce , 
		         	          const std::vector<Real>& Mass );


                void ThermostatIntegrateMLFF( std::vector<VecNSlice<Real> >& Pos ,
		                              std::vector<VecNSlice<Real> >& ForceDir ,
		                              std::vector<VecNSlice<Real> >& ForceCar ,
                                              std::vector<VecNSlice<Real> >& Velo ,
			                      Real dt ,
			                      const std::vector<std::vector<Real> >& Box , 
			                      const std::vector<std::vector<Real> >& RecBox , 
			                      MachineLearningForces GetForce , const std::vector<Real>& Mass );
		
                // particle insertion method alterning parameters
		void AddParticle( Real value );
		
		void ComputeLayerTemperature( const std::vector<VecNSlice<Real> >& Velos , 
				              const std::vector<Real>& Mass );
                
		void VelocityExchange( std::vector<VecNSlice<Real> >& Velos );






        private:
	        // declare field for Random number generator
   	        Ran Rand;  //initialization in constructor
		//Andersen 
		std::vector<Real> Sigma ;
		Real ColProb;       // collision property
		Real ColProbScaled; // collision probability times time step
		//~~~~~~~~~~~~~~~~~~

		// Nose-Hoover
		Real VirtMass ;     // Nose- Hoover virtual mass
		Real NoseChi  =  Real( 1.0 );
		Real NoseS    =  Real( 1.0 );
		Real NoseG;
		Real NoseError = Real( 1e-10 );
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~

		//Langevin
		Real LangevinFric;
		Real LangevinA;
		std::vector<Real> LangevinB;
		Real LangevinC;
		std::vector<VecNSlice<Real> > RandomForce;
		std::vector<Real> RandomForceData;
		//~~~~~~~~~~~~~
		
		// Bussi Donadio Parrinello
		Real BDPKinetic;   // desired kinetic energy
		Real BDPKineticAct;   // desired kinetic energy
		Real BDPDOF;       // degrees of freedom
		Real BDPDOFMOD2;       // degrees of freedom minus 1 / 2
		Real BDPFactor1;    // exponential factor for BDP
		Real BDPFactor2;    // exponential factor for BDP
		Real BDPFactor3;    // exponential factor for BDP
		Real BDPAlpha;
		bool evenDOF;      // odd or even computation of gamma random number
		// heat gradient BDP thermostat
		unsigned int EquiSteps;
		bool BDPHeatGrad ; // check if heat gradient is switched on or off
		std::vector<Real> BDPLayerAlpha; // alpha values for heat gradient computation
		std::vector<Real> BDPKineticLayer;  // desired kinetic energies in layers
		std::vector<Real> BDPKineticLayerAct; //actual kinetic energies in layers

		//Heat gradient variables
		HeatGradientStruc HeatVars;     // input parameters
		std::vector<std::vector<unsigned int> > HeatGradLink;  // link list assigning atoms to temperature
		std::vector<Real> TempGradient;  // temperature in regions of box
		
		std::vector<unsigned int> AMLinkList; // regions of box to thermostat



		
		
		


		unsigned int ThermoType;
		Real temperature;

		// leap frog algorigthm without thermostating
		// can be used in combination with the BDP thermostat
                void LeapFrogStepA( std::vector<VecNSlice<Real> >& Pos , 
		                    std::vector<VecNSlice<Real> >& Velo , 
		                    const std::vector<VecNSlice<Real> >& Force ,
			            const std::vector<Real>& Mass , Real dt );
                void LeapFrogStepB( std::vector<VecNSlice<Real> >& Velo ,
		                    const std::vector<VecNSlice<Real> >& Force ,
		                    const std::vector<Real>& Mass ,
				    Real dt ); 
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



		// velocity verlet algorithm with andersen thermostat
		void AndersenStepA( std::vector<VecNSlice<Real> >& Pos , 
				    const std::vector<VecNSlice<Real> >& Force ,
				    std::vector<VecNSlice<Real> >& Velo , Real dt , 
				    const std::vector<Real>& Mass );
		// call AndersenStepA then compute forces call AndersenStepB
		void AndersenStepB( std::vector<VecNSlice<Real> >& Velo ,
			            const std::vector<VecNSlice<Real> >& Force , 
		  		    Real dt , const std::vector<std::vector<Real> >& Box ,
				    const std::vector<Real>& Mass );
		// nose hoover thermostat (-> non stochastic) in velo-verlet leap frog representation
		void NoseHooverStepA( std::vector<VecNSlice<Real> >& Pos ,
				      const std::vector<VecNSlice<Real> >& Force ,
				      std::vector<VecNSlice<Real> >& Velo  ,
				      Real dt , const std::vector<std::vector<Real> >& Box , 
				      const std::vector<Real>& Mass );
		// propagate velocity from t+dt/2 to t+dt; compute forces inbetween step A and step B
		void NoseHooverStepB( const std::vector<VecNSlice<Real> >& Force ,
				      std::vector<VecNSlice<Real> >& Velo  ,
				      Real dt , const std::vector<std::vector<Real> >& Box , 
				      const std::vector<Real>& Mass );
		// Langevin dynamics with egeneralized langevin equation
		// first step in velocity verlet
		void LangevinStepA( std::vector<VecNSlice<Real> >& Pos ,
		                    const std::vector<VecNSlice<Real> >& Force ,
		                    std::vector<VecNSlice<Real> >& Velo ,
		                    Real dt ,
		                    const std::vector<std::vector<Real> >& Box ,
		                    const std::vector<Real>& Mass );
		// propagate the velocity from t+dt/2 to t + dt-> call after computing forces and step A
		void LangevinStepB( const std::vector<VecNSlice<Real> >& Force ,
		                    std::vector<VecNSlice<Real> >& Velo ,
		                    Real dt ,
		                    const std::vector<std::vector<Real> >& Box ,
		                    const std::vector<Real>& Mass );
		// generate random numbers mimicing the collisions or the drag in Langevin
	        void LangevinRandomNumbers( void );
		// Remove center of mass drift
		void RemoveCOMDrift( std::vector<VecNSlice<Real> >& Velo , const std::vector<Real>& Mass );


		//heat gradient helper functions
		void WriteHeatGradientLinkList( void );
                void MakeTemperatureGradient( void );
		void InitHeatGradientTemperatures( const std::vector<int>& NDirs );
		void InitHeatGradientLinkList( const std::vector<int>& NDirs );

		// bussi donadino parrinello thermostat
                void BDPComputeKineticEnergy( const std::vector<VecNSlice<Real> >& VeloCity ,
		                              const std::vector<std::vector<Real> >& Box , 
					      const std::vector<Real> & Mass );
                void BDPComputeKineticEnergyHeatGrad( 
		                        const std::vector<VecNSlice<Real> >& VeloCity ,
		                        const std::vector<std::vector<Real> >&Box , 
					const std::vector<Real> & Mass );

                void BDPRecomputeAlpha( void );
                void BDPRecomputeAlphaHeatGradAM( void );
                void BDPRecomputeAlphaHeatGrad( void );
		void BDPRescaleVelocities( std::vector<VecNSlice< Real> >& VeloCity );
                void BDPRescaleVelocitiesHeatGrad( std::vector<VecNSlice< Real> >& VeloCity );
                // rescale only some areas of box for active measure
                void BDPRescaleVelocitiesHeatGradAM( std::vector<VecNSlice< Real> >& VeloCity ); 
	

};
#endif
