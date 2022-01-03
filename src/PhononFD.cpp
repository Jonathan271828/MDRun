#include "PhononFD.h"



 std::vector<Real> extractVector( std::string data ){
     std::string Temp;
     std::vector<Real> RealData ;

     std::stringstream StrStream( trim( data.c_str() ) );
     while ( StrStream >> Temp ){
            RealData.push_back( ::atof( Temp.c_str() ) );
     }
     return RealData;
 }






 FDPhonons::FDPhonons( std::vector<std::vector<Real> > Box , std::vector<std::vector<Real> > RecBox ,
		       std::vector<int> NN , std::vector<Real> Dx ){

	 std::ifstream Qfile;
	 Qfile.open( "QPOINTS" );
	 std::string data;
	 delta = Dx;

	 unsigned int col = 0;
	 std::cout << "#######################################################\n" << std::endl;
	 std::cout << "Phonon computation by finite differences is swicthed on\n" << std::endl;
	 
	 for ( auto i = 0 ; i < Box.size(); i++ ){
		 //UnitCell.push_back( Box[i] / Real( NN[i] ) );
		 //RecUnitCell.push_back( RecBox[i] * Real( NN[i] ) );
		 UnitCell.push_back( Box[i] ); 
		 RecUnitCell.push_back( RecBox[i] );
	 }
	 std::cout << "Unit cell \n" << std::endl;
	 for ( auto i = 0 ; i < Box.size(); i++ ){
		 for ( auto j = 0 ; j < Box[i].size() ; j++ ){
			 printf( " %15.8f" , UnitCell[i][j] );
		 }
		 printf( "\n" );
	 }

	 std::cout <<  "Reciprocal unit cell \n" << std::endl;
	 for ( auto i = 0 ; i < Box.size(); i++ ){
		 for ( auto j = 0 ; j < Box[i].size() ; j++ ){
			 printf( " %15.8f" , RecUnitCell[i][j] );
		 }
		 printf( "\n" );
	 }

	 Real PI2 = M_PI * Real( 2.0 );

	 std::getline( Qfile , data );
	 std::getline( Qfile , data );
	 std::cout << "using QPOINTS in cartesian reciprocal coordinates times factor 2pi\n" << std::endl;
	 while ( std::getline( Qfile , data ) ){
		 Qpoints.push_back( extractVector( data ) );
		 //Qpoints[ col ]  =  PI2 * MatVec( RecUnitCell , Qpoints[ col ] );
		 Qpoints[ col ]  =  PI2 * Qpoints[ col ];
		 for ( auto i = 0 ; i < Qpoints[ col ].size() ; i++ ){
		     printf( "%15.8f" , Qpoints[col][i] / PI2 );
		 }
		     printf( "\n");
		 col++;
	 }

	 std::cout << "Finite difference is using displacements" << std::endl;
	 std::cout << "in direct coordinates of" << std::endl;
	 for ( auto i = 0 ; i < delta.size() ; i++ ){
		 std::cout << std::setw( 15 ) << std::setprecision( 8 ) << delta[i];
	 }
	 std::cout << std::endl;

	 std::cout << "#######################################################\n" << std::endl;
}


 // computing the force constant matrix by finite differences
 // first loop allocates the force constant matrix
 // then the second derivatives are computed for the pair interactions
 // since the matrix has to be symmetric only the upper half triagle is comuted and
 // the loewr trinagle is assigned according to
 void FDPhonons::ComputeForceConstantMatrix( std::vector<VecNSlice<Real> >& Pos , 
		        ComputeForces& GetForce , std::vector<std::vector<Real> >& Box , std::vector<Real>& Mass ){


	 // coded for single atom in unit cell
	 ForceConstMatrix.resize( 1 );
	 ForceConstMatrix[0].resize( 1 );
	 ForceConstMatrix[0][0].resize( Pos[0].size() );
	 for ( auto i = 0 ; i < Pos[0].size() ; i++ ){
		 ForceConstMatrix[0][0][i].resize( Pos[0].size() );
		 for ( auto j = 0 ; j < Pos[0].size(); j++ ){
			 ForceConstMatrix[ 0 ][ 0 ][ i ][ j ].resize( Pos.size() );
			 for ( auto k = 0 ; k < Pos.size() ; k++ ){
				 ForceConstMatrix[ 0 ][ 0 ][ i ][ j ][ k ].resize( Pos.size() );
			 }
		 }
	 }

	 DebugFile.open( "Debug_phon.out" );
	 std::ofstream TestFile;
	 TestFile.open( "Test.dat" );

	 unsigned int dim = Pos[0].size();

	 printf( "Computing force constant matrix in real space\n" );
	 std::clock_t t1 = std::clock();
	 for ( auto l = 0 ; l < Pos.size() ; l++ ){
		 for ( auto nn = 0 ; nn < GetForce.LinkList.NNList[l].size(); nn++ ){
			 unsigned int lp = GetForce.LinkList.NNList[l][nn];
		         for ( auto alpha = 0 ; alpha < dim ; alpha++ ){
				 for ( auto beta = 0 ; beta < dim ; beta++ ){
					 // off diagonal terms
					 ForceConstMatrix[ 0 ][ 0 ][ alpha ][ beta ][ l ][ lp ] =
					      	 Real( 2 ) * ComputeSecondDerivativeOD( Pos[l] , Pos[lp] , Box , 
					 			   alpha , beta , GetForce ) / 
					  	         sqrt( Mass[l] * Mass[ lp ] );
					 // l = lp elements
//					 ForceConstMatrix[ 0 ][ 0 ][ alpha ][ beta ][ l ][ l ] =
//					  	           ComputeSecondDerivativeDE( Pos , Box , alpha , 
//					 		   	    beta , GetForce , l ) / Mass[l];
//

					 //DebugFile << l << "   " << lp << "  " << alpha << "   ";
					 //for ( auto a = 0 ; a < dim ; a++ ){
					 //        DebugFile << Pos[l][a] << "   ";
					 //}
					 //for ( auto a = 0 ; a < dim ; a++ ){
					 //        DebugFile << Pos[lp][a] << "   ";
					 //}
					 //DebugFile << ForceConstMatrix[ 0 ][ 0 ][ alpha ][ beta ][ l ][ lp ]
					 //          << "   ";
					 //DebugFile << ForceConstMatrix[ 0 ][ 0 ][ alpha ][ beta ][ l ][ l ]
					 //          << std::endl;
				 }
			 }
		 }
		 DebugFile << std::endl;
	 }
	 std::clock_t t2 = std::clock();
	 printf( "Force constant computation took %16.4f sec\n" , Real( t2-t1 ) / CLOCKS_PER_SEC ) ;

	 //for ( auto j = 0 ; j < ForceConstMatrix[0][0][0][0].size() ; j++ ){ //l atom index
    	 //      for ( auto i = 0 ; i < ForceConstMatrix[ 0 ][ 0 ].size(); i++ ){ // xyz
	 //       	 for ( auto l = 0 ; l < ForceConstMatrix[0][0][0][0][0].size() ; l++  ){ // l prime
	 //       	         for ( auto k = 0 ; k < ForceConstMatrix[0][0][0].size() ; k++ ){ //x'y'z'
	 //       			 DebugFile << ForceConstMatrix[0][0][i][k][j][l] << "   ";
	 //       		 }
	 //       	 }
	 //        DebugFile << std::endl;
	 //        }
	 //}
	 //throw;
 }


  // helper function to compute second derivative of
  // used potential potential function has to be of the form ( PosA , PosB , Box )
  // only works for off diagonal elements
  Real FDPhonons::ComputeSecondDerivativeOD( VecNSlice<Real> PosA , VecNSlice<Real> PosB , 
		                           std::vector<std::vector<Real>> Box , unsigned int alpha,
					   unsigned int beta , ComputeForces GetForce ){

	  Real InitA = PosA[ alpha ];
	  Real InitB = PosB[ beta ];

	  // sift A + delta and B + delta   factor 0.5 because of forward backward difference
	  PosA[ alpha ] = InitA + Real( 0.5 ) * delta[ alpha ];
	  PosB[ beta ]  = InitB + Real( 0.5 ) * delta[ beta ];
	  Real Epp = GetForce.SupplyMorsePairOnlyPB( PosA , PosB , Box );
	  // sift A - delta and B - delta
	  PosA[ alpha ] = InitA - Real( 0.5 ) * delta[ alpha ];
	  PosB[ beta ]  = InitB - Real( 0.5 ) * delta[ beta ];
	  Real Emm = GetForce.SupplyMorsePairOnlyPB( PosA , PosB , Box );
	  // sift A - delta and B + delta
	  PosA[ alpha ] = InitA - Real( 0.5 ) * delta[ alpha ];
	  PosB[ beta ]  = InitB + Real( 0.5 ) * delta[ beta ];
	  Real Emp = GetForce.SupplyMorsePairOnlyPB( PosA , PosB , Box );
	  // sift A + delta and B - delta
	  PosA[ alpha ] = InitA + Real( 0.5 ) * delta[ alpha ];
	  PosB[ beta ]  = InitB - Real( 0.5 ) * delta[ beta ];
	  Real Epm = GetForce.SupplyMorsePairOnlyPB( PosA , PosB , Box );
	  PosA[alpha] =  InitA;
	  PosB[beta]  =  InitB;
	  GetForce.SupplyMorsePairOnlyPB( PosA , PosB , Box );
	  return ( Epp + Emm - Epm - Emp ) / delta[ alpha ] / delta[ beta ];
  }


  // compute second derivative for diagonal element
  Real FDPhonons::ComputeSecondDerivativeDE( std::vector<VecNSlice<Real> >& Pos , 
		                             std::vector<std::vector<Real> >& Box,
		                             unsigned int alpha , unsigned int beta , ComputeForces& GetForce , 
	       	                             unsigned int AtNr ){
	  Real InitA = Pos[ AtNr ][ alpha ];
	  Real InitB = Pos[ AtNr ][ beta ];
	  if ( alpha != beta ){
		  // alpha + h , beta + h
		  Pos[AtNr][ alpha ] = InitA + Real( 0.5 ) * delta[ alpha ];
		  Pos[AtNr][ beta ] = InitB  + Real( 0.5 ) * delta[ beta ];
		  Real Epp = GetForce.SupplyMorseNNEOnlyPB( AtNr , Pos , Box );
		  // alpha - h , beta - h
		  Pos[AtNr][ alpha ] = InitA - Real( 0.5 ) * delta[alpha];
		  Pos[AtNr][ beta ] = InitB  - Real( 0.5 ) * delta[ beta ];
		  Real Emm = GetForce.SupplyMorseNNEOnlyPB( AtNr , Pos , Box );
		  // alpha - h , beta + h
		  Pos[AtNr][ alpha ] = InitA - Real( 0.5 ) * delta[alpha];
		  Pos[AtNr][ beta ]  = InitB  + Real( 0.5 ) * delta[ beta ];
		  Real Emp = GetForce.SupplyMorseNNEOnlyPB( AtNr , Pos , Box );
		  // alpha + h , beta - h
		  Pos[AtNr][ alpha ] =  InitA + Real( 0.5 ) * delta[alpha];
		  Pos[AtNr][ beta ]  =  InitB  - Real( 0.5 ) * delta[ beta ];
		  Real Epm = GetForce.SupplyMorseNNEOnlyPB( AtNr , Pos , Box );

		  Pos[ AtNr ][ alpha ] = InitA ;
		  Pos[ AtNr ][ beta ] = InitB ;

		  return ( Epp + Emm - Epm - Emp ) / delta[alpha] / delta[ beta ];
	  }
	  else{
		  Real E0 =  GetForce.SupplyMorseNNEOnlyPB( AtNr , Pos , Box );
		  Pos[ AtNr ][ alpha ] = InitA + 0.5*delta[ alpha ]; 
		  Real Ep = GetForce.SupplyMorseNNEOnlyPB( AtNr , Pos , Box );
		  Pos[ AtNr ][ alpha ] = InitA - 0.5*delta[ alpha ];
		  Real Em = GetForce.SupplyMorseNNEOnlyPB( AtNr , Pos , Box );
		  Pos[ AtNr ][ alpha ] =  InitA;
		  Pos[ AtNr ][ beta ]  =  InitB;
		  return ( Ep + Em - 2.0 * E0 ) / delta[ alpha ] / delta[ alpha ];
	  }
  }






 void FDPhonons::ComputeDynamicMatrix( std::vector<VecNSlice<Real> > Pos , std::vector<std::vector<Real> > Box ){


	 DynamicMatrix.resize( Qpoints.size() );
	 for ( auto i = 0 ; i < Qpoints.size() ; i++ ){
		 DynamicMatrix[i].resize( ForceConstMatrix.size()*ForceConstMatrix[0][0].size() ,
				 ForceConstMatrix.size()*ForceConstMatrix[0][0].size() ); // initialized with zeros
	 }


	 printf( "Computing dynamical matrix by \n" );
	 printf( "fourier transforming the force constant matrix\n" );

	 std::clock_t t1 = std::clock();

	 // vector filled with zeros to center atomic positions
	 // around the origin
	 std::vector<Real> Center;
	 Center.resize( Pos[0].size() );

	 std::ofstream factorfile;
	 factorfile.open( "TestFile.dat" );
	 //
	 //std::ofstream Test23;
	 //Test23.open( "Test23.dat" );


	 // normalization constant for the 
	 // fourier transform
	 Real norm = Real( 1.0 ) ;//Real( Pos.size() );
	 unsigned int dim = Pos[0].size();
	 for ( auto qq = 0 ; qq < Qpoints.size() ; qq++ ){
		 for ( auto k = 0 ; k < ForceConstMatrix.size(); k++ ){ // dummy loop
			 for ( auto kp = 0 ; kp < ForceConstMatrix.size(); kp++ ){ // dummy loop
				 for ( auto alpha = 0 ; alpha < dim ; alpha++ ){//xyz
					 unsigned int A  =  k*dim  +  alpha;
					 for ( auto beta = 0 ; beta < dim ; beta++ ){
						 unsigned int B = kp*dim  +  beta;
						 for ( auto l = 0 ; l < Pos.size() ; l++ ){
							 for ( auto lp = 0 ; lp < Pos.size() ; lp++ ){
								 std::vector<Real> TP =
									 get_nearest_image( Pos[l] , Pos[lp] );
								 TP = TP - Pos[l];
								// std::cout << TP[0] << "   "
								//	   << TP[1] << "   "
								//	   << TP[2] << std::endl;
								 CReal factorB =  
								        exp( CReal( Real( 0 ) , DotProduct( TP , Qpoints[qq] ) ) );  
							 
								// factorfile << factorB.real() << "   " 
							 	//    << factorB.imag() << std::endl;



								 DynamicMatrix[ qq ]( A , B )  =
								 DynamicMatrix[ qq ]( A , B ) + 
								 ForceConstMatrix[k][kp][alpha][beta][l][lp] * factorB;
							 }
						 }
					         DynamicMatrix[ qq ]( A , B )  =  DynamicMatrix[ qq ]( A , B ) / 
							                   //CReal( sqrt( Real( Pos.size() ) ) , 0 );
							                   CReal( Real( Pos.size() ) , 0 );
					 }
				 }
			 }
		 }
         }
	 std::clock_t t2 = std::clock();
	 printf( "Fourier transformation took %16.4f sec\n" , Real( t2-t1 ) / CLOCKS_PER_SEC );
 }


 // computing eigenvalues and eigenvectors of dynamic matrix
 // stored in cpp vector DynamicMatrix, every enetry of vector is
 // contain the elements of dynamic matrix for certain q-vector
 // in form of Eigen Matrix
 void FDPhonons::ComputeEigenVectorsValues( void ){


	 // check for the eigensolver -> it gives the correct result
	 // example can be found here : https://eigen.tuxfamily.org/dox/classEigen_1_1ComplexEigenSolver.html
	 /*****************************************************
	 Eigen::Matrix4cd test;
	 test << CReal(-0.211,0.68) ,  CReal(0.108,-0.444) ,  CReal(0.435,0.271)  , CReal(-0.198,-0.687),
		 CReal(0.597,0.566) ,  CReal(0.258,-0.0452),  CReal(0.214,-0.717) , CReal(-0.782,-0.74),
		 CReal(-0.605,0.823),  CReal(0.0268,-0.27) ,  CReal(-0.514,-0.967), CReal(-0.563,0.998),
		 CReal(0.536,-0.33) ,  CReal(0.832,0.904)  ,  CReal(0.608,-0.726) , CReal(0.678,0.0259);
	 ComplexSolver.compute( test );
	 std::cout << ComplexSolver.eigenvalues() << std::endl;
	 ******************************************************/


	 EigenVectors.resize( Qpoints.size() );
	 EigenValues.resize( Qpoints.size() );
	 for ( auto i = 0 ; i < DynamicMatrix.size() ; i++ ){
		 ComplexSolver.compute( DynamicMatrix[i] , Eigen::ComputeEigenvectors );
		 EigenValues[i].resize( ComplexSolver.eigenvalues().size() );
		 EigenVectors[i].resize( ComplexSolver.eigenvalues().size() );
		 for ( auto j = 0 ; j < EigenValues[i].size() ; j++ ){
			 EigenValues[i][j] = ComplexSolver.eigenvalues()(j);
			 for ( auto k = 0 ; k < ComplexSolver.eigenvectors().rows();k++ ){
			     EigenVectors[i][j].push_back( ComplexSolver.eigenvectors()(j,k) );
			 }
		 }
	 }
 }

 /*
  *
  * Write real and imaginary phonon dispersion
  * curves to files along the supplied q-path
  *
  */

 void FDPhonons::WriteOutput( void ){

	 std::vector<Real> Qpath;
	 Qpath.resize( Qpoints.size() );
	 Qpath[0] = Real( 0.0000000 );
	 for ( auto i = 1 ; i < Qpoints.size() ; i++ ){
		 Qpath[ i ] =  VecNorm( Qpoints[ i ] - Qpoints[ i - 1 ] ) + Qpath[ i - 1 ]; 
	 }

	 std::ofstream DispFile;
	 std::ofstream DispFileImag;
	 DispFile.open( "PhononDispReal.out" );
	 DispFileImag.open( "PhononDispImag.out" );
	 for ( auto i = 0 ; i < Qpoints.size() ; i++ ){
		 DispFile << std::setprecision( 8 ) << std::setw( 15 ) << Qpath[i] ;
		 DispFileImag << std::setprecision( 8 ) << std::setw( 15 ) << Qpath[i] ;
		 for ( auto j = 0 ; j < EigenValues[i].size(); j++ ){
		     DispFile << std::setprecision( 8 ) << std::setw( 15 ) << EigenValues[i][j].real() ;
		     DispFileImag << std::setprecision( 8 ) << std::setw( 15 ) << EigenValues[i][j].imag() ;
		 }
		 DispFile << std::endl;
		 DispFileImag << std::endl;
	 }
	 DispFile.close();
	 DispFileImag.close();

	 // write dynmat file

 }


 /*
  * Routine writes output in the quantum espresso format such that
  * it can be used for the phq code as projectors
  * on the harmonic phonon eigenmodes
  *
  */
 void FDPhonons::WriteQuantumEspressoOutput( void ){

	 std::ofstream Outfile;
	 Outfile.open( "dyn.out" );

	 const unsigned int widthQ = 14;
	 const unsigned int precQ = 9;
	 const unsigned int widthVec = 10;
	 const unsigned int precVec = 6;


	 Outfile.setf( std::ios::fixed , std::ios::floatfield );

	 for ( auto i = 0 ; i < EigenValues.size(); i++ ){
	         Outfile.precision( precQ );
		 Outfile << std::endl;
		 Outfile <<"      Dynamical  Matrix in cartesian axes\n";
		 Outfile << std::endl;
		 Outfile << "     q = ( ";
		 for ( auto j = 0 ; j < Qpoints[i].size(); j++ ){
			 Outfile << std::setw( widthQ ) << Qpoints[i][j];
		 }
		 Outfile << " )" << std::endl;
		 Outfile << std::endl;
		 Outfile << " ";
		 for ( auto i = 0 ; i < 74 ; i++ ){
			 Outfile << "*";
	         }
		 Outfile << std::endl;
		 Outfile.precision( precVec );
		 for ( auto j = 0 ; j < EigenValues[i].size(); j++ ){
		         Outfile << "      freq (    " << std::setw( 1 ) << j+1 <<") =       ";
			 Outfile << std::setw( 8 ) << EigenValues[i][j].real() * eVtTHz
			       	 << " [THz] =       ";
			 Outfile << std::setw( 8 ) << EigenValues[i][j].real() * eVtIcm 
			         << " [cm-1]" << std::endl;
			 Outfile << " (";
			 for ( auto k = 0 ; k < EigenVectors[i][j].size(); k++ ){
			     Outfile << std::setw( widthVec ) << EigenVectors[i][j][k].real() 
				     << std::setw( widthVec ) << EigenVectors[i][j][k].imag();
		         }
			 Outfile << " )" << std::endl;
		 }
		 for ( auto i = 0 ; i < 74 ; i++ ){
			 Outfile << "*";
		 }
		 Outfile << std::endl;
		 Outfile << std::endl;
		 Outfile.flush();
	 }
	 Outfile.close();
 }



 /*
  * Routine writes output in a very simple format
  * first line contains number of qpoints, and number of branches per point
  * empty line
  * eigenvalue 1
  * eigenvector 1
  * eigenvalue 2
  * eigenvector 2
  * . . .
  * . . .
  * . . .
  * empty line next qpoint
  * next qpoint
  * it can be used for the eval_rouinte
  * for harmonic eigenmode projection phonon eigenmodes
  *
  */
 void FDPhonons::WriteSimpleOutput( void ){

	 std::ofstream Outfile;
	 Outfile.open( "dyn_eval.out" );

	 const unsigned int widthQ = 14;
	 const unsigned int precQ = 9;
	 const unsigned int widthVec = 10;
	 const unsigned int precVec = 6;
	 Real PI2 = M_PI * Real( 2.0 );


	 Outfile.setf( std::ios::fixed , std::ios::floatfield );

         Outfile << Qpoints.size() << "   " << EigenValues[0].size() << "  " << std::endl;;
	 for ( auto i = 0 ; i < EigenValues.size(); i++ ){
	         Outfile.precision( precQ );
		 for ( auto j = 0 ; j < Qpoints[i].size(); j++ ){
			 Outfile << std::setw( widthQ ) << Qpoints[i][j] / PI2;
		 }
		 Outfile << std::endl;
		 Outfile.precision( precVec );
		 for ( auto j = 0 ; j < EigenValues[i].size(); j++ ){
			 //Outfile << std::setw( 8 ) << EigenValues[i][j].real() << std::endl;
			 for ( auto k = 0 ; k < EigenVectors[i][j].size(); k++ ){
			     Outfile << std::setw( widthVec ) << EigenVectors[i][j][k].real();
				     //<< std::setw( widthVec ) << EigenVectors[i][j][k].imag();
		         }
			 Outfile << std::endl;
			 Outfile << std::endl;
		 }
		 Outfile << std::endl;
		 Outfile.flush();
	 }
	 Outfile.close();
 }



 // restore the acoustic sum rule for the dynamic lattice
 void FDPhonons::AcousticSumRule( unsigned int dim ){


	 Eigen::MatrixXcd Ksum;
	 Ksum.resize( dim , dim );

	 unsigned int N = dim -1 ;

	 CReal NATS = CReal( DynamicMatrix[0].rows() ) / CReal( 3 ) ;

	 for ( auto qq = 0 ; qq < DynamicMatrix.size() ; qq++ ){
		 if ( VecNorm( Qpoints[ qq ] ) < 1e-5 ){
			 for ( auto i = 0 ; i < DynamicMatrix[qq].rows(); i=i+3 ){
				 Ksum.setZero();
				 for ( auto j = 0 ; j < DynamicMatrix[qq].cols(); j=j+3 ){
					 Ksum = Ksum + DynamicMatrix[qq].block( i, j, dim, dim );
				 }
				 Ksum = Ksum / NATS;
				 for ( auto j = 0 ; j < DynamicMatrix[qq].cols(); j=j+3 ){
					 DynamicMatrix[qq].block(i,j,dim,dim) = 
						 DynamicMatrix[qq].block(i,j,dim,dim) - Ksum;
				 }
			 }
		 }
	 }
 }


 // impose symmetry on the dynamic matrix
 void FDPhonons::Symmetrize( unsigned int dim ){
	 for ( auto qq = 0 ; qq < DynamicMatrix.size(); qq++ ){
		 for ( auto i = 0 ; i < DynamicMatrix[qq].rows() ; i++ ){
			 for ( auto j = 0 ; j < DynamicMatrix[qq].cols(); j++ ){
				 CReal Value = ( DynamicMatrix[qq]( i , j ) + 
						 DynamicMatrix[ qq ]( j , i ) ) * Real( 0.5 );
				 DynamicMatrix[ qq ]( i , j ) = Value;
				 DynamicMatrix[ qq ]( j , i ) = std::conj( Value );
			 }
	         }
	 }
 }

void FDPhonons::PrintDynMatrix( unsigned int QN ){

	 for ( auto i = 0 ; i < DynamicMatrix[ QN ].rows() ; i++ ){
		 for ( auto j = 0 ; j < DynamicMatrix[ QN ].cols() ; j++ ){
			 DebugFile << " " << DynamicMatrix[ QN ]( i , j ).real() << "  "
				   << DynamicMatrix[ QN ]( i , j ).imag() << "  ";
		 }
		 DebugFile << std::endl;
	 }
}


void FDPhonons::WriteForceMatrix( void ){


	std::vector<std::ofstream> Outfiles;
	Outfiles.resize( 3 );
	Outfiles[0].open( "XX.dat" );
	Outfiles[1].open( "YY.dat" );
	Outfiles[2].open( "ZZ.dat" );


	 for ( auto k = 0 ; k < ForceConstMatrix[0][0][0].size() ; k++ ){ //x'y'z'
		 for ( auto j = 0 ; j < ForceConstMatrix[0][0][0][0].size() ; j++ ){ //l atom index
			 for ( auto l = 0 ; l < ForceConstMatrix[0][0][0][0][0].size() ; l++  ){ // l prime
	        			 Outfiles[ k ] << ForceConstMatrix[0][0][k][k][j][l] << "   ";
	        	 }
	                 Outfiles[ k ] << std::endl;
	         }
	         Outfiles[ k ].close();
	 }
}



 void  FDPhonons::main( std::vector<VecNSlice<Real> >& Pos , ComputeForces& GetForce , 
		        std::vector<std::vector<Real> >& Box , std::vector<Real>& Mass ){

	 ComputeForceConstantMatrix( Pos , GetForce , Box , Mass );
         WriteForceMatrix();
	 ComputeDynamicMatrix( Pos , Box );
	 PrintDynMatrix( 1 );
	 //AcousticSumRule( Pos[0].size() );
	 //Symmetrize( Pos[0].size() );
	 //PrintDynMatrix( 4 );
	 ComputeEigenVectorsValues();
	 WriteOutput();
	 WriteQuantumEspressoOutput();
	 WriteSimpleOutput();
	 printf( "Output finite differences written \n" );
 }









