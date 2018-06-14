#include "lib_iris.hpp"

#define NB_IMG_MAX 250
int main ( int argc, char ** argv ) 
{
	c_learning obj1;
	
	
	if ( argc > 1 )
	{
		if ( 	!	strcmp( argv[1], "--help" ) 	||
				!	strcmp( argv[1], "-H" ) )
		{
			obj1.help(cout, 1);
			cout << "Author: Valérian Némesin." << endl;
			return 0;
		}
	}
	
	if ( obj1.setup( argc, argv ) )
		return 1;
		
		
	if ( obj1.segment( ) )
		return 1;
		
	return 0;
	
	
}
