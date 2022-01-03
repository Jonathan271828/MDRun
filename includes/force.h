#include <omp.h>

#include <stdio.h>
#ifndef _PeriodicBound
#include "PeriodicBoundary.h"
#endif
#ifndef _VectorTools
#include "MatMulVec.h"
#endif
#ifndef _LinkList
#include "LinkList.h"
#endif

#ifndef _MDSPEC
#include "atoms.h"
#endif


#ifdef USE_DOUBLES
typedef double Real;
#else
typedef float Real;
#endif




#ifndef _Forces
#define _Forces


class OnSitePotential{

	private:
		std::vector<Real> Coupling;
		std::vector<VecNSlice<Real> > EquiPos;
		std::vector<Real> EquiPosArray;


		unsigned int Factorial( unsigned int );

        public:
		OnSitePotential( void ){};
	        void OnSitePotentialInit( const Real CP , const int Order , 
			const std::vector<VecNSlice<Real> > & EquiPos );
		Real ComputeEnergyOnSite( const std::vector<VecNSlice<Real> >&Pos ,
				          const std::vector<std::vector<Real> >& Box );
		Real EnergySingleAtom( const Real norm );
		void ComputeForceOnSite( const std::vector<VecNSlice<Real> > & Pos , 
				         const std::vector<std::vector<Real> >& Box ,
					 const std::vector<std::vector<Real> >& RecBox ,
				         std::vector<VecNSlice<Real> >& Force );
                std::vector<Real> ComputeForceOnSiteSingle( 
		                       const VecNSlice<Real> & Pos ,
		                       const std::vector<std::vector<Real> >& Box ,
				       const std::vector<std::vector<Real> >& RecBox , 
				       const int index );
                Real ComputeEnergySingleAtomOnSite( const VecNSlice<Real>& Pos ,
		                                    const std::vector<std::vector<Real> >&Box,
	     	                                    const int index );
};


 /*
  *
  *
  * imposes an on site Mexican force to the
  * atoms; by setting the second parameter zero
  * one can use this class to impose harmonic force
  *
  */

class DoubleWellPotential{

	public:
		Real DoubleWellParA;         //parameters for potential second order
		Real DoubleWellParB;         //parameters for potential 4th order
		Real DoubleWellParCCart;
		Real DoubleWellEquiDistCart;  //equilibrium distance for double well potential
		std::vector<std::vector<Real> > CoP;  // center of potential
		Real Normfactor;
	        Real ExpFactor; 
	        Real HarmFactor;

		DoubleWellPotential( void ){};

		// setting parameters for the double well potential
                void SetDoubleWellParameters( Real A , Real B , Real C , 
		                              Real EquiDistCart );
		// compute pair force taking into account periodic boundary conditions
                std::vector<Real> DoubleWellForcePairPB( const VecNSlice<Real> & PosA , 
	  	                                         const VecNSlice<Real> & PosB ,
							 const std::vector<std::vector<Real>>& Box,
							 const std::vector<std::vector<Real>>& RecBox );
		// compute pair force without periodic boundary conditions
                std::vector<Real> DoubleWellForcePair( const VecNSlice<Real> & PosA , 
		                                       const VecNSlice<Real> & PosB ,
		     			               const std::vector<std::vector<Real>>& Box,
						       const std::vector<std::vector<Real>>& RecBox );
		// compute potential energy between a pair of atoms PosA PosB with taking
		// into account periodic boundary conditions
                Real DoubleWellEpotPairPB( const VecNSlice<Real>& PosA , const VecNSlice<Real>& PosB ,
 		                           const std::vector<std::vector<Real> >& Box );
		// compute potential energy between a pair of atoms PosA and PosB taking
		// NO periodic boundary conditions
                Real DoubleWellEpotPair( const VecNSlice<Real>& PosA , const VecNSlice<Real>& PosB ,
		                         const std::vector<std::vector<Real> >& Box );
};




/*
 *
 *
 * class object gives the user the possibility
 * to compute harmonic forces between atoms
 *
 *
 */



class HarmonicForce{


	public:
		Real HarmonicPotPar ;   // coupling constant
	        Real HarmonicEquiDistCar;  // distance between harmonicly coupled atoms	cartesian
		std::vector<Real> HarmonicEquiDists; // equilibrium distances between atoms


                HarmonicForce( void ){};
		~HarmonicForce( void ){};
		Real ComputePotential( const std::vector<VecNSlice<Real> >& Pos ,
                                       const std::vector<std::vector<int> >& LinkList ,
				       const std::vector<std::vector<Real> >& Box );

		void SetHarmonicForceParameters( Real x , Real EquiDistCar );

		// compute pair forces with individual parameters
                std::vector<Real> HarmonicForcePair( const VecNSlice<Real>& PosA , const VecNSlice<Real>& PosB ,
	        			             const std::vector<std::vector<Real>> & Box ,
						     const std::vector<std::vector<Real>> & RecBox );
                std::vector<Real> HarmonicForcePairPB( const VecNSlice<Real>& PosA ,
		                                       const VecNSlice<Real>& PosB , 
						       const std::vector<std::vector<Real>> & Box ,
						       const std::vector<std::vector<Real>> & RecBox );


		// compute potential energies with varying parameters
                Real HarmonicEpotPair( std::vector<Real> PosA , std::vector<Real> PosB ,
		                       std::vector<std::vector<Real>> Box );
                Real HarmonicEpotPairPB( const VecNSlice<Real>& PosA , const VecNSlice<Real>& PosB ,
		                         const std::vector<std::vector<Real> >& Box );
};




 //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


 /*
  *
  *
  * Morse potential between atoms
  *
  *
  */


 class MorsePotential{

	 public:
		 Real MorseD;   // depth of morse potential at equilibrium
		 Real MorsealphaCart ;
		 Real MorseEquiDistCart;   // equilibrium distance cartesian

		 MorsePotential( void ){};
		 void SetMorseForceParameters( Real x , Real descentCart , Real EquiDistCar );


		 // force functions
                 std::vector<Real> MorseForcePairPB( const VecNSlice<Real>& PosA , const VecNSlice<Real>& PosB ,
		                                     const std::vector<std::vector<Real> >& Box , 
		                                     const std::vector<std::vector<Real> >& RecBox );

                 std::vector<Real> MorseForcePair( const VecNSlice<Real>& PosA , const VecNSlice<Real>& PosB, 
		                                   const std::vector<std::vector<Real> >& Box , 
		                                   const std::vector<std::vector<Real> >& RecBox );


		 // potential functions
		 Real MorseEpotPair( const VecNSlice<Real>& PosA , const VecNSlice<Real>& PosB ,
                             const std::vector<std::vector<Real> >& Box );
		 Real MorseEpotPairPB( const VecNSlice<Real>& PosA , const VecNSlice<Real>& PosB ,
                               const std::vector<std::vector<Real> >& Box );

 };













//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


 /*
  *
  * force object that inherits from all the other
  * force classes. Module combines all forces
  *
  */


 class ComputeForces: public HarmonicForce,DoubleWellPotential,
	                     MorsePotential,OnSitePotential{

	 public:
		 // fixed list implementation
		 ComputeForces( const std::vector<VecNSlice<Real> >& Pos , unsigned int NN , 
	  	                const std::vector<std::vector<Real> >& Box );
		
		 // verlet list update force implementation
		 ComputeForces( const std::vector<VecNSlice<Real> >& PosCar ,
                                const std::vector<VecNSlice<Real> >& PosDir ,
				Real cutoff , 
		 	        const std::vector<std::vector<Real> >& Box );
				        

		 ComputeLinkLists LinkList;


		 void SetUpOnSitePotential( const Real CP , const int Order ,
				const std::vector<VecNSlice<Real> >& PosEqui );

		 // init harmonic force/potential
		 void SetUpHarmonicPart( Real x , Real EquiDistCar );
		 // init Morse force / potential
		 void SetUpMorsePart( Real D , Real alphaCart , Real EquiDistCar );
		 // init double well pair potential for disordered system
                 void SetUpDoubleWellPart( Real A , Real B , Real C , 
		                           Real EquiDist );

		 // double well potential forces
                 void SupplyDoubleWellOnly( const std::vector<VecNSlice<Real> >& Pos ,
		                           const std::vector<std::vector<Real> > & Box ,
		                           const std::vector<std::vector<Real> > & RecBox ,
   		                           std::vector<VecNSlice<Real> >& force );
		 // harmonic forces
		 //void SupplyHarmonicOnly( const std::vector<VecNSlice<Real> >& Pos ,
		 //	                  std::vector<VecNSlice<Real> >& force );
                 void SupplyHarmonicOnly( const std::vector<VecNSlice<Real> >& Pos ,
                                          const std::vector<std::vector<Real> >& Box ,
                                          const std::vector<std::vector<Real> >& RecBox ,
   		                          std::vector<VecNSlice<Real> >& force );
		 // Morse force
                 void SupplyMorseOnly( const std::vector<VecNSlice<Real> >& Pos , 
		                       const std::vector<std::vector<Real> > & Box ,
		                       const std::vector<std::vector<Real> > & RecBox ,
		                       std::vector<VecNSlice<Real> >& force );

		 // compute force array for on site potential
		 void SupplyOnSiteForceOnly( const std::vector<VecNSlice<Real> >& Pos ,
           				     const std::vector<std::vector<Real> >& Box ,     
                                             const std::vector<std::vector<Real> > &RecBox , 
				             std::vector<VecNSlice<Real> >& force );

		 // compute force arising due to harmonic nearest neighbour potential and
		 // on site potential
		 void SupplyHarmonicAndOnSiteForce( const std::vector<VecNSlice<Real> >& Pos ,
						    const std::vector<std::vector<Real> >& Box ,
						    const std::vector<std::vector<Real> > &RecBox ,
						    std::vector<VecNSlice<Real> >& force );


		 // harmonic potential
		 Real PotentialEnergyHarmonic( const std::vector<VecNSlice<Real> >& Pos ,
				               const std::vector<std::vector<Real> >& Box );

		 // double well potential
		 Real PotentialEnergyDoubleWell( const std::vector<VecNSlice<Real> >& Pos ,
		  		                 const std::vector<std::vector<Real> >& Box );
		 // Morse potential
		 Real PotentialEnergyMorse( const std::vector<VecNSlice<Real> >& Pos ,
				            const std::vector<std::vector<Real> >& Box );
		 // compute Morse potential for certain atomic pair
                 Real SupplyMorsePairOnlyPB( const VecNSlice<Real>& PosA , const VecNSlice<Real>& PosB , 
		                             const std::vector<std::vector<Real> >& Box );
		 // computes the energy for atom with index AtNr in the morse potential
                 Real SupplyMorseNNEOnlyPB( unsigned int AtNr , const std::vector<VecNSlice<Real> >& Pos , 
		                            const std::vector<std::vector<Real> >& Box );

		 // supply on site energy as Nth order polynomial
		 Real SupplyOnSiteEnergyOnly( const std::vector<VecNSlice<Real> >& Pos ,
				              const std::vector<std::vector<Real> >& Box );
 
		  // Harmonic nearest neighbour coupling plus double well potential coupling
		  Real PotentialEnergyHarmonicAndDoubleWell( 
				  const std::vector<VecNSlice<Real> >& Pos , 
		                  const std::vector<std::vector<Real> >& Box );

		 // mixed force functions
		 std::vector<std::vector<Real> > HarmonicMexican( std::vector<std::vector<Real>> Pos );

		 // mixed force harmonic and mexican
		 void SupplyHarmonicAndOnSite( const std::vector<VecNSlice<Real> >& Pos ,
				                     std::vector<VecNSlice<Real> >& force );
		 // harmonic potential energy and on site potential
		 Real PotentialEnergyHarmonicAndOnSite( const std::vector<VecNSlice<Real> >& Pos , 
				                        const std::vector<std::vector<Real> >& Box );

		 // morse potential for a single atom
		 Real SupplyMorseSingleAtom( const VecNSlice<Real>& PosNew , 
				             const std::vector<VecNSlice<Real> >& Pos , 
                                             const std::vector<unsigned int>& List , 
					     const std::vector<std::vector<Real> >& Box );
                 std::vector<Real> SupplyMorseForceSingleAtom( const VecNSlice<Real>& PosNew ,
		                                             const std::vector<VecNSlice<Real> >& Pos , 
		                                             const std::vector<unsigned int>& List ,
							     const std::vector<std::vector<Real>>& Box ,
							     const std::vector<std::vector<Real>>& RecBox );

		 // filling morse potential energy array
		 void ComputeMorsePotentialEnergyArray( const std::vector<VecNSlice<Real> >& Pos ,
                                                        const std::vector<std::vector<Real> >& lattice ,  
                                                        std::vector<Real>& Epot );
	
 };
#endif
