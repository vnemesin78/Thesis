#include "lib_iris.hpp"
#include <cstdio>

int main( int argc, char ** argv )
{
	if ( argc == 1 )
	{
		cout << "Error : missing argument(s)" << endl;
		return 1;
	}

	if ( 	! strcmp( argv[1], "--help" )	||
			! strcmp( argv[1], "-H" ) )
	{
		cout << "HELP" << endl;
		cout << " Arg[1] : Src Rep " << endl;
		cout << " Arg[2] : Save Rep" << endl;
		cout << " Arg[3 - n] : Parameter files" << endl;
		cout << "Author: Valérian Némesin." << endl;
		return 0; 
	}	
	
	if ( argc < 3 )
	{
		cout << "Error : missing argument(s)" << endl;
		return 1;
	}
	
	api_parameters params;
	for ( unsigned int i = 3; i < argc; ++ i )
		params.load( argv[i] );
	
	c_fusion_iris_template toto;
	
	if ( toto.setup( params ) )
		return 1;
	
	if ( ! toto.fusion( argv[1] ) )
		toto.save( argv[2] );
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
}
