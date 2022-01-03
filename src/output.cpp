

#ifndef _OUTPUT
#include "output.h"
#endif





 
 void ErrorMessage( std::string text , std::string routine ){

	 std::cout << "###########################" << std::endl;
	 std::cout << "~~~~~~~~~~ ERROR ~~~~~~~~~~" << std::endl;
	 std::cout << "In routine" << routine << std::endl;
	 std::cout << text                    << std::endl;
	 std::cout << "~~~~~~~~~~ ERROR ~~~~~~~~~~" << std::endl;
	 std::cout << "###########################" << std::endl;

 }




 void WriteOutput::XdatInit( std::vector<std::vector<Real> > lattice , std::vector<std::string> Spec ,
		 std::vector <int> Natoms ){

	 fprintf( file , "Molecular dynamics output\n" );
	 fprintf( file , "   1   \n" );



	 for ( auto j = 0 ; j < lattice.size(); j++ ){
		 for ( auto i = 0 ; i < lattice[j].size(); i++ ){
			 fprintf( file , "   %f" , lattice[j][i] );
		 }
	         fprintf( file , " \n" );
	 }
	 for ( auto i = 0 ; i < Spec.size() ; i++ ){
		 fprintf( file , "  %s" , Spec[i].c_str() );
         }
	 fprintf( file , "\n" );
	 for ( auto i = 0 ; i < Natoms.size() ; i++ ){
		 fprintf( file , "  %i" , Natoms[i] );
         }
	 fprintf( file , "\n" );
	 Counter = 1;
 }

 void WriteOutput::AddStruc( std::vector<VecNSlice<Real> > Pos ){

	 fprintf( file , "Direct    =      %i\n" , Counter );


	 for ( auto i = 0 ; i < Pos.size() ; i++ ){
		 for ( auto j = 0 ; j < Pos[i].size(); j++ ){
			 fprintf( file , "   %15.8f" , Pos[i][j] );
	         }
		 fprintf( file , "\n" );
	 }
	 Counter = Counter + 1;
 }


 void WriteOutput::LammpsWrite( std::vector<VecNSlice<Real> > Pos , std::vector<std::vector<Real> > lattice ,
		               std::vector<int> Spec , unsigned int Natoms ){
	 fprintf( file , "ITEM: TIMESTEP\n" );
	 fprintf( file , "0\n");
	 fprintf( file , "ITEM: NUMBER OF ATOMS\n" );
	 fprintf( file , "%i \n" , Pos.size() );
	 fprintf( file , "ITEM: BOX BOUNDS pp pp pp\n" );

	 for ( auto i = 0 ; i < lattice.size() ; i++ ){
		 fprintf( file , "   %15.8f   %15.8f  \n" , Real( 0 ) , lattice[i][i] );
	 }
	 fprintf( file , "ITEM: ATOMS id type xs ys zs\n" );
	 for ( auto i = 0 ; i < Pos.size() ; i++ ){
		 fprintf( file , "  %i  %i" , i+1 , Spec[i] );
		 //std::vector<Real> temp = MatVec( lattice , Pos[i] );
		 for ( auto j = 0 ; j < Pos[i].size() ; j++ ){
		         //fprintf( file , "   %15.8f   " , temp[j] );
		         fprintf( file , "   %15.8f   " , Pos[i][j] );
		 }
		 fprintf( file , "\n");
     }
         //fflush( file );     // slowes the code but therefore output is written immidiately
 }

 void WriteOutput::PhQWriteStruc( std::vector<VecNSlice<Real> > Pos , bool first , unsigned int Number ){

	 if ( not first ){
		 fprintf( file , "\nmd_step  =  %i\n" , Number );
	 }
	 else{
		 fprintf( file , "total_step  =  %i\n" , Number );
		 fprintf( file , "\n" );
         }
	 fprintf( file , "atomic_positions\n" );
	 for ( auto i = 0 ; i < Pos.size(); i++ ){
	         fprintf( file , "Cu" );
		 for ( auto j = 0 ; j < Pos[i].size() ; j++ ){
			 fprintf( file , "%18.10f" , Pos[i][j] );
		 }
		 fprintf( file , "\n" );
	 }
         //fflush( file );     // slowes the code but therefore output is written immidiately
 }




 void WriteOutput::WriteVector( std::vector<Real> data ){

	 fprintf( file , "   %i" , Counter );
	 for ( auto i = 0 ; i < data.size(); i++ ){
		 fprintf( file , "   %f" , data[i] );
	 }
	 fprintf( file , "\n" );
	 fflush( file );     // slowes the code but therefore output is written immidiately
	 Counter = Counter + 1;
 }

