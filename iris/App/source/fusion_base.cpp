#include "lib_iris.hpp"
#include <dirent.h>
unsigned int get_nb_classes ( const char * rep )
{
	unsigned int i = 0;
	DIR * dir;
	dir = opendir( rep );
	if ( !dir )
	{
		cout << "Can't open " << rep << "!" << endl;
		return 1;
	}
	struct dirent * ent;
	while ( ( ent = readdir( dir ) ) != NULL )
	{
		if ( 	strcmp( ent->d_name, ".." ) 			&& 
				strcmp( ent->d_name, "." )	 			)
		{
			stringstream oss;
			oss << rep << "/" << ent->d_name;
			if ( is_dir( oss.str().c_str() ) )
				++ i;
		}
		
	}
	closedir(dir);
	return i;
}

int get_class_names ( 	const char * rep,
						string ** names,
						unsigned int & nb_classes )
{
	nb_classes = get_nb_classes ( rep ); 
	if ( !nb_classes )
	{
		cout << "Error : No sub directories" << endl;
		return 1;
	}
	
	*names = new string[nb_classes];
	
	unsigned int i = 0;
	DIR * dir;
	dir = opendir( rep );
	if ( !dir )
	{
		cout << "Can't open " << rep << "!" << endl;
		return 1;
	}
	struct dirent * ent;
	while ( ( ent = readdir( dir ) ) != NULL )
	{
		if ( 	strcmp( ent->d_name, ".." ) 			&& 
				strcmp( ent->d_name, "." )	 			)
		{
			stringstream oss;
			oss << rep << "/" << ent->d_name;
			if ( is_dir( oss.str().c_str() ) )
			{
				(*names)[i] = ent->d_name;
				++ i;
			}
			
		}
		
	}
	closedir(dir);	
	return 0;
}


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
	for ( int i = 3; i < argc; ++ i )
		params.load( argv[i] );
	
	c_fusion_iris_template toto;
	if ( toto.setup( params ) )
		return 1;
	
	
	unsigned int nb_classes;
	string * names;
	
	if ( get_class_names ( 	argv[1],
							&names,
							nb_classes ) )
		return 1;

#ifdef __linux
	mkdir ( argv[2], 014777 );
#elif __APPLE__ 
	mkdir ( argv[2], 014777 );
#elif _WIN32
	_mkdir( argv[2] );
#elif _WIN64
	_mkdir( argv[2] );
#else
	#error
#endif


	for ( unsigned int i = 0; i < nb_classes; ++ i )
	{
		stringstream oss1, oss2;
		oss1 << argv[1] << "/" << names[i];
		oss2 << argv[2] << "/" << names[i];
		cerr << oss1.str() << endl;
		toto.fusion( oss1.str().c_str() );
		toto.save( oss2.str().c_str()  );
		int q = system("clear");
		if ( q == 1 )
			cout << "Ubuntu!" << endl;
		cout << "(" << ( ( 1000 * i ) / nb_classes ) / 10.0 << "/100)" << endl;
	}
	
	return 0;
	
	
}
