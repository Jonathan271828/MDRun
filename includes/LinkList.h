#include <stdio.h>
#include <vector>
#include <limits>
#include <math.h>
#include <omp.h>


#ifndef _VectorTools
#include "MatMulVec.h"
#endif
#ifndef _PeriodicBound
#include "PeriodicBoundary.h"
#endif


#ifndef _ANALYZE
#include "analyze.h"
#endif

#ifndef _MDSPEC
#include "atoms.h"
#endif





#ifdef USE_DOUBLES
typedef double Real;
#else
typedef float Real;
#endif





#ifndef _LinkList
#define _LinkList
 class ComputeLinkLists
 {

	private:
		Real minVal = std::numeric_limits<Real>::min();
		std::vector<std::vector<unsigned int> > NNBoxList;
		std::vector<std::vector<unsigned int> > VerletNNList;
		std::vector<std::vector<Real> > VerletNNDists;
		std::vector<std::vector<Real> > VerletNNTravel;
		Real VerletCutoff;
		Real VerletDist;
		// helper routine for ComputeNNListPeriodicBound to sort the indices 
		// update the nearest neighbour field for every atom properly
		// during list computation
		void UpdateField( std::vector<std::vector<Real> > & data , Real Dist , unsigned int NewNum, 
				       unsigned int Index ); 
                void MakeBoxes( std::vector<std::vector<Real> > Box , Real Cutoff , 
				        std::vector<std::vector<Real> > InverseBox );
		std::vector<std::vector<std::vector<Real> > > BoxBounds;
		std::vector<std::vector<unsigned int> > BoxList;
		std::vector<std::vector<unsigned int> > BoxToAtom;
	        std::vector<unsigned int> NBoxes;
                void DetermineNeighbouringBoxes( void );
                unsigned int ReturnBoxNumber( const unsigned int x , 
				              const unsigned int y , const unsigned int z );
                std::vector<unsigned int> GetNeighbourIndexSingleBox( unsigned int x , 
				                                      unsigned int y , 
								      unsigned int z );
                unsigned int ReturnBoxPos1D( const int x , const unsigned int Nmax , const int Add );
                bool CheckVerletUpdate( const std::vector<VecNSlice<Real> >& Velos , 
				        Real timestep , Real Cutoff );
                void UpdateList( Real Cutoff );
                void RezeroVerletNNTravel( unsigned int N1 , unsigned int N2 );
                void UpdateVerletNNDists( std::vector<VecNSlice<Real> > Pos , 
				          const std::vector<std::vector<Real> > Box );

		void ComputeAverageNN( void );



	public: 
		std::vector<std::vector<int> > NNList;
		ComputeLinkLists( void ){};
		// computes the first NN nearest neighbours of an atom
		// based on a simple sorting scheme
		void ComputeNNListPeriodicBound( const std::vector<VecNSlice<Real> >& Pos , unsigned int N , 
				                         const std::vector<std::vector<Real> >& Box );

                // N^2 scaling algorithm
		void ComputeNNListPeriodicBound( const std::vector<VecNSlice<Real> >& Pos , 
		                                 const std::vector<std::vector<Real> >& Box , Real cutoff );
                void PrintNNList( void );
                void ComputeNNListPeriodicBoundInit( const std::vector<VecNSlice<Real> >& PosCar , 
	        			             const std::vector<VecNSlice<Real> >& PosDir , 
						     Real Cutoff ,
                                                     const std::vector<std::vector<Real> >& Box ); 
                                                     
                void AssignParticlesToBoxes( const std::vector<VecNSlice<Real> >& Pos , Real Cutoff , 
	        			     const std::vector<std::vector<Real> >& Box );
                void ComputeNearestNeighborListBox( const std::vector<VecNSlice<Real> >& Pos ,
	        		                    const std::vector<std::vector<Real> >& Box ,
						    Real Cutoff );
   	        			 
                void UpdateNNVerletList( const std::vector<VecNSlice<Real> >& PosDir ,
				         const std::vector<VecNSlice<Real> >& PosCar, 
	        			 const std::vector<VecNSlice<Real> >& Velos , 
	        	                 const Real timestep , const Real Cutoff , 
	        			 const std::vector<std::vector<Real> >& Box );
                void ComputeNearestNeighborN2( const std::vector<std::vector<Real> >& Pos , 
                                               const std::vector<std::vector<Real> >& Box , Real Cutoff );
 };
#endif
