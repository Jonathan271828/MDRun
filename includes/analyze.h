


#include <vector>
#include <stdio.h>
#include <iostream>
#include <math.h>
#include <string>
#include <omp.h>
#include <algorithm>


#ifndef _MDSPEC
#include "atoms.h"
#endif
#ifndef __VectorTools
#include "MatMulVec.h"
#endif
#ifndef _OUTPUT
#include "output.h"
#endif
#ifndef _PeriodicBound
#include "PeriodicBoundary.h"
#endif
#ifndef _INTEGRATE
#include "SimpsIntegration.h"
#endif

#ifndef _MDSPEC
#include "atoms.h" 
#endif





#ifdef USE_DOUBLES
typedef double Real;
#else
typedef float Real;
#endif



#ifndef _ANALYZE
#define _ANALYZE

 unsigned int computeBin( Real x , Real dx , unsigned int Max );



 Real compute_spherical_theta( std::vector<Real> x );

 Real compute_spherical_phi( std::vector<Real> x );


 void SphericalSinCorrection( std::vector<std::vector<Real> > & data , int dim = 0 );



 Real ComputeKineticEnergy( std::vector<VecNSlice<Real> > Velo , std::vector<Real> Mass , std::vector<Real>& Ekin );


 Real ComputeTemperature( Real Ek , int Nf );


 Real ComputePressure( std::vector<VecNSlice<Real> > Pos , std::vector<VecNSlice<Real> > force ,
	             const std::vector<std::vector<Real> > Box , const Real T ,
	             const std::vector<std::vector<int> > List );

#pragma omp declare reduction( vec_Real_plus : std::vector<Real> : \
                               std::transform( omp_out.begin(), omp_out.end(), omp_in.begin(), omp_out.begin(), std::plus<Real>())) \
                               initializer(omp_priv = decltype(omp_orig)(omp_orig.size()))



 class PairDistFunction
     {
       private:
	     std::vector<Real> distribution;
	     std::vector<Real> xaxis ;
	     Real dx ;
	     unsigned int Nbins;

	     std::vector<Real> PairDistribution;
	     std::vector<Real> PairAxis;
	     Real PairDx;
	     unsigned int PairBins;



       public:
	     std::vector<std::vector<Real>> Center;

		 PairDistFunction( void ){};
		 void DistCenter( Real xmax , std::vector<VecNSlice<Real> > Cent , int Nmax = 100 );
		 void PairDistInit( Real xmax , int Nmax = 100 );
		 void ComputeDistToCenter( std::vector<VecNSlice<Real>> Pos ,
                                           std::vector<std::vector<Real>> lattice );
                 void ComputePairDistribution( std::vector<VecNSlice<Real> > Pos , 
			                       std::vector<std::vector<Real> > Box );
		 void WriteOutputCOM( std::string fname );
		 void WriteOutputPD( std::string fname );



 };


 class AngularDist3d{

	 private:
		 unsigned int Ntheta;
		 unsigned int Nphi;
		 Real dtheta;
		 Real dphi;
		 std::vector<std::vector<Real>> distribution;

	public:
		 AngularDist3d( unsigned int N1 , unsigned int N2 );
                 void ComputeAngularCOMCentered( std::vector<VecNSlice<Real> > Pos ,
		                                 std::vector<std::vector<Real> > COM ,
		       			         std::vector<std::vector<Real> > lattice );
		 void WriteOutput( std::string fname );


 };

class AnalyzePotential{



	private:
		std::vector<VecNSlice<Real> > EquiPos;
		std::vector<Real> EquiPosArray;
		std::vector<std::vector<Real> > Histograms;
		std::vector<std::vector<unsigned int> > Counter;
		std::vector<std::vector<Real> > ForceHistogram;
		unsigned int Nbins = 500;
		Real min  =  Real( 0 );
		Real max  =  Real( 5 );
		Real dx   =  ( max - min ) / Real( Nbins );

		std::vector<VecNSlice<Real> > Center;
		std::vector<Real> CenterArray;
		std::vector<VecNSlice<Real> > Delta;
	        std::vector<Real> DeltaArray;

		std::vector<Real> AxisArray;
		std::vector<VecNSlice<Real>> Axis;
		std::vector<Real> FaceArray;
		std::vector<VecNSlice<Real>> Face;
		std::vector<Real> RoomArray;
		std::vector<VecNSlice<Real>> Room;


                void AssignVector( std::vector<Real>& data , std::vector<VecNSlice<Real> >& Pointers , unsigned int dim );
		void MakeAxisDirections( void );
		void MakeFaceDirections( void );
		void MakeRoomDirections( void );
		unsigned int CheckWhichDirection( const std::vector<Real>& delta );

	public:
		AnalyzePotential( const std::vector<VecNSlice<Real> >& Pos );
		void PotentialAnaMain( const std::vector<VecNSlice<Real> >& Pos , const std::vector<Real>& Epot,
				       const std::vector<std::vector<Real> >& lattice ,
				       const std::vector<VecNSlice<Real> >& Force );
		void PotentialAnaFinalize( void );

};

#endif
