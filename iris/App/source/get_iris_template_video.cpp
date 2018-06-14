#include "lib_iris.hpp"

int main( int argc, char ** argv )
{
	int end = 0;
	if ( argc == 1 )
	{
		cout << "Error : missing argument(s)" << endl;
		return 1;
	}

	if ( 	! strcmp( argv[1], "--help" )	||
			! strcmp( argv[1], "-H" ) )
	{
		cout << "HELP" << endl;
		cout << " Arg[1] : Image filename " << endl;
		cout << " Arg[2] : Save directory " << endl;
		cout << " Arg[3 - n] : Parameter files" << endl;
		cout << "Author: Valérian Némesin." << endl;
		return 0; 
	}	
	
	if ( argc < 3 )
	{
		cout << "Error : missing argument(s)" << endl;
		return 1;
	}
	
	c_get_iris_template_pthread obj;

	if ( obj.setup( argc - 3, argv + 3, argv[2], &cout, true) )
	{
		cout << "prob1" << endl;
		return 1;
	}

	
	if ( ! obj.run( argv[1], end ) )
		obj.save( argv[2] );	
	else
		cout << "prob" << endl;
	return 0;
}
