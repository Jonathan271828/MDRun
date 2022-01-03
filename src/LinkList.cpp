#ifndef _LinkList
#include "LinkList.h"
#endif



 // helper routine for ComputeNNListPeriodicBound 
 // to update the tempData field properly
 void ComputeLinkLists::UpdateField( std::vector<std::vector<Real> > & data , Real Dist ,
		                       unsigned int NewNum , unsigned int Index ){

	 std::vector<std::vector<Real> > temp;
	 temp.resize( data.size() - Index );
	 // store old values
	 for ( auto i = 0 ; i < temp.size() ;i++ ){
		 temp[i].resize( data[i+Index].size() );
		 temp[i][0] = data[i+Index][0];
		 temp[i][1] = data[i+Index][1];
	 }
	 data[Index][0] = Dist;
	 data[Index][1] = Real( NewNum );
	 for ( auto i = 0 ; i < temp.size()-1; i++ ){
		 data[Index+i+1][0] = temp[i][0];
		 data[Index+i+1][1] = temp[i][1];
     }
 }




// compute ordinary link list with and assign to every pos the first N nearest neighbor positions
// in direct coordinates
  void ComputeLinkLists::ComputeNNListPeriodicBound( const std::vector<VecNSlice<Real> >& Pos , unsigned int N , 
		                                     const std::vector<std::vector<Real> >& Box ){


	  std::vector<std::vector<Real> >  tempData;
	  tempData.resize( N );
	  for ( auto i = 0 ; i < N ; i++ ){
		  tempData[i].resize( 2 );
		  tempData[i][0]  =  Real( 1e200 );
		  tempData[i][1]  =  Real( 0 );
	  }

	  printf( "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" );
	  printf( "Nearest neighbour connection table \n" );
	  NNList.resize( Pos.size() );
	  for ( auto i = 0 ; i < Pos.size() ; i++ ){
		  NNList[ i ].resize( N );
		  for ( auto j = 0 ; j < Pos.size() ; j++ ){
		      if ( i == j ) continue;
		      std::vector<Real> delta = get_nearest_image( Pos[i] , Pos[j] ) - Pos[i];
		      Real norm = VecNorm( MatVec( Box , delta ) );
		      for ( auto k = 0 ; k < N ; k++ ){
			      if ( norm <  tempData[k][0] ){ 
				      UpdateField( tempData , norm , j , k );
				      break;
			      }
		      }
	      }

		  for ( auto k = 0 ; k < N ; k++ ){
		      NNList[i][k]  =  int( tempData[k][1] );
		      printf( " %i " , NNList[i][k] );
		      tempData[k][0]  =  Real( 1e200 );
		  }
		  printf( "  \n");
	  }
	  //PrintNNList();
	  printf( "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" );
  }
	

  /// compute nearest neighbor list with a N^2 algorithm
  void ComputeLinkLists::ComputeNNListPeriodicBound( const std::vector<VecNSlice<Real> >& Pos , 
				                     const std::vector<std::vector<Real> >& Box ,
						     Real cutoff ){

	  NNList.resize( Pos.size() );
	  for ( auto i = 0 ; i < Pos.size(); i++ ){
		  unsigned int NN = 0;
		  NNList[i].resize( 0 );
		  for ( auto j = 0 ; j < Pos.size() ; j++ ){
			  if ( i == j ) continue;
		      std::vector<Real> delta = get_nearest_image( Pos[i] , Pos[j] ) - Pos[i];
			  Real norm = VecNorm( MatVec( Box , delta ) );
			  if ( norm < cutoff ){
				  NNList[i].push_back( j );
				  NN += 1; 
			  }
		  }
      }
  }


  // Write nearest neighbor list to screen
  void ComputeLinkLists::PrintNNList( void ){

	  printf( "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" );
	  printf( "Nearest neighbour connection table \n" );
	  for ( auto i = 0 ; i < NNList.size(); i++ ){
		  for ( auto j = 0 ; j < NNList[i].size() ; j++ ){
			  printf( " %i " , NNList[ i ][ j ] );
          }
		  printf( "  \n");
	  }
	  printf( "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" );
  }



  // compute link list in cell decomposition of simulation box
  void ComputeLinkLists::ComputeNNListPeriodicBoundInit( const std::vector<VecNSlice<Real> >& PosCar ,
		                                         const std::vector<VecNSlice<Real> >& PosDir ,
							 Real Cutoff ,
                                                         const std::vector<std::vector<Real> >& Box ){
	  // make this in constructor
	  VerletCutoff =  Real( 1.5 ) * Cutoff;
	  VerletDist   =  VerletCutoff - Cutoff;
	  AssignParticlesToBoxes( PosCar , Cutoff , Box );
         //////////////////

	  DetermineNeighbouringBoxes(); 
	  RezeroVerletNNTravel( PosCar.size() , PosCar[0].size() );
	  ComputeNearestNeighborListBox( PosDir , Box , VerletCutoff );
	  UpdateList( Cutoff );
  }



  // replace this bullshit by a binning routine this is much more efficient
  void ComputeLinkLists::AssignParticlesToBoxes( const std::vector<VecNSlice<Real> >& Pos , Real Cutoff , 
		                                         const std::vector<std::vector<Real> >& Box ){

	  NBoxes.resize( Box.size() );
	  for ( auto i = 0 ; i < Box.size() ; i++ ){
		  NBoxes[i]  = ( unsigned int )( VecNorm( Box[i] ) / Cutoff );
	  }

	  NNBoxList.resize( Product( NBoxes ) );
	  BoxToAtom.resize( NNBoxList.size() );
	  BoxList.resize( Pos.size() );
	  // resizing boxes
	  for ( auto i = 0 ; i < BoxToAtom.size() ; i++ ){
		  BoxToAtom[i].resize( 0 );
	  }


	  for ( auto i = 0 ; i < Pos.size() ; i++ ){
		  BoxList[i].resize( Pos[i].size() );
		  // atoms -> box 
		  for ( auto j = 0 ; j < Pos[i].size() ; j++ ){
		      BoxList[ i ][ j ] = computeBin( Pos[ i ][ j ] , Cutoff , NBoxes[ j ] );
		  }
		  // box -> atom
		  unsigned int NBox  =  ReturnBoxNumber( BoxList[ i ][ 0 ] , 
				                         BoxList[ i ][ 1 ] , 
							 BoxList[ i ][ 2 ] );
	 	  BoxToAtom[ NBox ].push_back( i );
	  }
  }



  // defines the conncetion table for neighboring boxes
  void ComputeLinkLists::DetermineNeighbouringBoxes( void ){

	  for ( auto i = 0 ; i < NBoxes[0] ; i++ ){
		  for ( auto j = 0 ; j < NBoxes[1] ; j++ ){
			  for ( auto k = 0 ; k < NBoxes[2] ; k++ ){
			      unsigned int col = ReturnBoxNumber( i , j , k );
				  NNBoxList[ col ] = GetNeighbourIndexSingleBox( i , j , k );
			  }
		  }
	  }
  }

  // return the position of a box when neighboring cell x+add is taken
  unsigned int ComputeLinkLists::ReturnBoxPos1D( const int x , const unsigned int Nmax , const int Add ){

	  int Result = x + Add;
	  if ( Result < 0 ){
		  Result = Nmax - 1;
	  }else if ( Result >= Nmax ) {
		  Result = 0;
	  }
	  return Result;
  }


  // returns indices for neighbouring boxes of box x,y,z
  std::vector<unsigned int> ComputeLinkLists::GetNeighbourIndexSingleBox( unsigned int x , 
		                                                          unsigned int y , 
									  unsigned int z ){

	  std::vector<unsigned int> List;
	  for ( auto i = -1 ; i < 2 ; i++ ){
		  unsigned int Nx  =  ReturnBoxPos1D( x , NBoxes[0] , i );
		  for ( auto j = -1 ; j < 2 ; j++ ){
		      unsigned int Ny  =  ReturnBoxPos1D( y , NBoxes[1] , j );
			  for ( auto k = -1 ; k < 2 ; k++ ){
		          unsigned int Nz  =  ReturnBoxPos1D( z , NBoxes[2] , k );
				  List.push_back( ReturnBoxNumber( Nx , Ny , Nz ) );
			  }
		  }
	  }
	  return List;
  }



  // give integer positions of box and return box Index
  unsigned int ComputeLinkLists::ReturnBoxNumber( const unsigned int x , 
		                                  const unsigned int y , 
						  const unsigned int z ){
	  return z + y * NBoxes[ 2 ] + x * ( NBoxes[ 1 ] * NBoxes[ 2 ] );
  }



  // routine needs direct coordinates on input and computes the VerletNNList
  void ComputeLinkLists::ComputeNearestNeighborListBox( const std::vector<VecNSlice<Real> >& Pos , 
		                                        const std::vector<std::vector<Real> >& Box ,
							Real Cutoff ){


	  VerletNNList.resize( Pos.size() );
	  VerletNNDists.resize( Pos.size() );

	  for ( auto i = 0 ; i < Pos.size() ; i++ ){
		  VerletNNList[ i ].resize( 0 );
		  VerletNNDists[ i ].resize( 0 );
		  unsigned int BoxIndx = ReturnBoxNumber( BoxList[ i ][ 0 ] ,
				                          BoxList[ i ][ 1 ] , 
							  BoxList[ i ][ 2 ] );
		  // loop over neighbouring boxes and box itself 
		  for ( auto  j = 0 ; j < NNBoxList[ BoxIndx ].size() ; j++ ){
			  // loop ovver atoms in neighbor boxes
			  unsigned int NBoxIdx =  NNBoxList[ BoxIndx ][ j ];
			  for ( auto k = 0 ; k < BoxToAtom[ NBoxIdx ].size() ; k++ ){
			      Real norm = VecNorm( MatVec( Box , get_nearest_image( 
							   Pos[i] , 
							   Pos[ BoxToAtom[ NBoxIdx ][ k ] ] ) - Pos[i] ) );
			      if ( norm < Cutoff && i != BoxToAtom[ NBoxIdx ][ k ] ) {
					  //std::cout << norm << std::endl;
				      VerletNNList[ i ].push_back( BoxToAtom[ NBoxIdx ][ k ] );
				      VerletNNDists[i].push_back( norm );
			      }
			  }
		  }
	  }
  }
 


  void ComputeLinkLists::ComputeNearestNeighborN2( const std::vector<std::vector<Real> >& Pos , 
		                                   const std::vector<std::vector<Real> >& Box ,
						   Real Cutoff ){

	  NNList.resize( Pos.size() );
	  for ( auto i = 0 ; i < Pos.size() ; i++ ){
		  NNList[i].resize( 0 );
		  unsigned int BoxIndx = ReturnBoxNumber( BoxList[ i ][ 0 ] ,
				                          BoxList[ i ][ 1 ] , 
							  BoxList[ i ][ 2 ] );
		  // loop over atoms in same box itself
		  for ( auto j = 0 ; j < BoxToAtom[ BoxIndx ].size() ; j++ ){
		  	  Real norm = VecNorm( MatVec( Box , 
					       get_nearest_image( Pos[i] , 
					       Pos[ BoxToAtom[ BoxIndx ][ j ] ] ) - Pos[i] ) );
			  if ( norm < Cutoff && i != BoxToAtom[ BoxIndx ][ j ] ) {
				  NNList[ i ].push_back( BoxToAtom[ BoxIndx ][ j ] );
			  }
		  }

		  // loop over neighbouring boxes 
		  for ( auto  j = 0 ; j < NNBoxList[ BoxIndx ].size() ; j++ ){
			  // loop ovver atoms in neighbor boxes
			  unsigned int NBoxIdx =  NNBoxList[ BoxIndx ][ j ];
			  for ( auto k = 0 ; k < BoxToAtom[ NBoxIdx ].size() ; k++ ){
			      Real norm = VecNorm( MatVec( Box , 
						   get_nearest_image( Pos[i] , 
						   Pos[ BoxToAtom[ NBoxIdx ][ k ] ] ) - Pos[i] ) );
			      if ( norm < Cutoff && i != BoxToAtom[ NBoxIdx ][ k ] ) {
					  NNList[ i ].push_back( BoxToAtom[ NBoxIdx ][ k ] );
				  }
			  }
		  }
	  }
	  //PrintNNList();
  }



  //
  //  checks if atoms moved trough verlet shell
  bool ComputeLinkLists::CheckVerletUpdate( const std::vector<VecNSlice<Real> >& Velos ,
		                            Real timestep , Real Cutoff ){

	  
	  for ( auto i = 0 ; i < VerletNNTravel.size(); i++ ){
		  VerletNNTravel[ i ] = VerletNNTravel[i] + Velos[i];
		  Real norm = VecNorm( VerletNNTravel[i] * timestep );
		  if ( norm > VerletDist ){
			  return true;
		  }
	  }
	  return false;
  }




  void ComputeLinkLists::UpdateList( Real Cutoff ){

	  NNList.resize( VerletNNList.size() );
	  for ( auto i = 0 ; i < VerletNNList.size() ; i++ ){
	      NNList[i].resize( 0 );
	      for ( auto j = 0 ; j < VerletNNList[i].size() ; j++ ){
			  if ( VerletNNDists[i][j] < Cutoff ){
				  NNList[i].push_back( VerletNNList[i][j] );
			  }
		  }
	  }
  }



  void ComputeLinkLists::UpdateVerletNNDists( std::vector<VecNSlice<Real> > Pos ,
		                              const std::vector<std::vector<Real> > Box ){

	  for ( auto i = 0 ; i < VerletNNList.size(); i++ ){
		  for ( auto j = 0 ; j < VerletNNList[i].size() ; j++ ){
			  VerletNNDists[i][j] = VecNorm( MatVec( Box , 
						get_nearest_image( Pos[ i ] , 
						Pos[ VerletNNList[ i ][ j ] ] ) - Pos[i] ) );
		  }
	  }
  }


  void ComputeLinkLists::UpdateNNVerletList( const std::vector<VecNSlice<Real> >& PosDir , 
                                             const std::vector<VecNSlice<Real> >& PosCar ,
		                             const std::vector<VecNSlice<Real> >& Velos ,
					     const Real timestep , 
					     const Real Cutoff , 
					     const std::vector<std::vector<Real> >& Box ){


	  if ( !CheckVerletUpdate( Velos , timestep , VerletDist ) && PosDir.size() == VerletNNDists.size() ){
	      UpdateVerletNNDists( PosDir , Box );
	  }else{
	    std::cout << "Updating Verlet-List" << std::endl;
	    AssignParticlesToBoxes( PosCar , Cutoff , Box );
	    ComputeNearestNeighborListBox( PosDir , Box , VerletCutoff );
	    RezeroVerletNNTravel( PosDir.size() , PosDir[0].size() );
	    //ComputeNNListPeriodicBound(PosDir , Box , Cutoff );
	  }
	  UpdateList( Cutoff );
  }



  void ComputeLinkLists::RezeroVerletNNTravel( unsigned int N1 , unsigned int N2 ){
	  
      VerletNNTravel.resize( N1 );
//#pragma parallel for
	  for ( auto i = 0 ; i < N1 ; i++ ){
	      VerletNNTravel[i].assign( N2 , Real( 0 ) );
	  }
  }



 void ComputeLinkLists::ComputeAverageNN( void ){

	 unsigned int NN = 0;

	 for ( auto i = 0 ; i < NNList.size() ; i++ ){
		 NN = NN + NNList[i].size();
	 }
	 std::cout << "NNN  " << float( NN ) / float( NNList.size() ) << std::endl;
 }





  /**************************************************
   **************************************************
   *************  JUNK SECTION **********************
   **************************************************
   **************************************************/

  void ComputeLinkLists::MakeBoxes( std::vector<std::vector<Real> > Box , Real Cutoff ,
		                            std::vector<std::vector<Real> > InverseBox ){
	  
	  
	  std::vector<unsigned int> NBoxes;
	  NBoxes.resize( Box.size() );

	  for ( auto i = 0 ; i < Box.size() ; i++ ){
		  NBoxes[i]  = ( unsigned int )( VecNorm( Box[i] ) / Cutoff );
	  }

	  BoxBounds.resize( Product( NBoxes ) );
	  // compute box centers
	  Real CutD2 = Cutoff / Real( 2 ); 
	  unsigned int NN = 0;
	  for ( auto i = 0 ; i < NBoxes[0] ; i++ ){
		  for ( auto j = 0 ; j < NBoxes[1] ; j++ ){
			  for ( auto k = 0 ; k < NBoxes[2] ; k++ ){
				  BoxBounds[NN].resize( NBoxes.size() );
				  for ( auto l = 0 ; l < NBoxes.size() ; l++ ){
					  BoxBounds[NN][l].resize( 2 );
				  }

				  std::vector<Real> r(3);
				  r[0] =  Cutoff * Real( i ) + CutD2;
				  r[1] =  Cutoff * Real( j ) + CutD2;
				  r[2] =  Cutoff * Real( k ) + CutD2;
				  r = MatVec( InverseBox , r );
				  BoxBounds[ NN ][ 0 ][ 0 ]  =  r[0];
				  BoxBounds[ NN ][ 1 ][ 0 ]  =  r[1];
				  BoxBounds[ NN ][ 2 ][ 0 ]  =  r[2];
				  
				  r[0] =  Cutoff * Real( i ) - CutD2;
				  r[1] =  Cutoff * Real( j ) - CutD2;
				  r[2] =  Cutoff * Real( k ) - CutD2;
				  r = MatVec( InverseBox , r );
				  BoxBounds[ NN ][ 0 ][ 1 ]  =  r[0];
				  BoxBounds[ NN ][ 1 ][ 1 ]  =  r[1];
				  BoxBounds[ NN ][ 2 ][ 1 ]  =  r[2];
				  NN += 1;
			  }
		  }
	  }
  }



