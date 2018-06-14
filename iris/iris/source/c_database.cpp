#include "c_database.hpp"
#include "distances.hpp"
#include "utils.hpp"
c_database :: c_database ( void )
{
	initialize();
}

c_database :: c_database ( 	const char * rep,
							const char * iris_code_file,
							const char * fragility_map_file,
							ostream * stream )
{
	initialize();
	setup ( rep, 
			iris_code_file,
			fragility_map_file,
			stream );
}
int c_database :: setup ( 	const char * rep,
							const char * iris_code_file,
							const char * fragility_map_file,
							ostream * stream )
{
	free();
	initialize();
	err_stream = stream;
	
	//Test de la validité des arguments d'entrée
	if ( 	!iris_code_file 	|| 
			!fragility_map_file || 
			!rep 				)
	{
		if ( err_stream )
			*err_stream << "Error : Assertion if ( !rep || !iris_code_file || !fragility_map_file ) failed!" << endl;
		return 1;
	}

	//Chargement de la base de données
	if ( load( 	rep, 
				iris_code_file, 
				fragility_map_file ) )
	{
		if ( err_stream )
			*err_stream << "Error : Failed to load iris database in " << rep << "!" << endl;
		return 1;	
	}

	return 0;
}

c_database :: ~c_database()
{
	free();
	initialize();
}

const IplImage * c_database :: iris_code( unsigned int id ) const
{
	if ( id < _nb_classes )
		return _iris_codes[id];
	else
		return NULL;
	
}

unsigned long long * c_database :: iris_code_bis( unsigned int id) const
{
	if ( id < _nb_classes )
		return _iris_codes_bis[id];
	else
		return NULL;
	
}
			
			
unsigned long long * c_database :: mask( unsigned int id) const
{
	if ( id < _nb_classes )
		return _masks[id];
	else
		return NULL;
}
			
unsigned long long * c_database :: fragile_bits( unsigned int id) const
{
	if ( id < _nb_classes )
		return _fragile_bits[id];
	else
		return NULL;
}


const IplImage * c_database :: fragility_map( unsigned int id ) const
{
	if ( id < _nb_classes )
		return _fragility_maps[id];
	else
		return NULL;
	
}


const char * c_database :: name( unsigned int id ) const
{
	if ( id < _nb_classes )
		return _names[id].c_str();
	else
		return NULL;
}

const char * c_database :: class_name( unsigned int id ) const
{
	if ( id < _nb_classes )
		return _class_names[id].c_str();
	else
		return NULL;
}

void c_database :: free()
{
	if ( _names )
		delete[] _names;
	if ( _class_names )
		delete[] _class_names;
	if ( _iris_codes )
	{
		for ( unsigned int i = 0; i < _nb_classes; ++ i )
		{
			cvReleaseImage ( _iris_codes + i );
		}
		delete[] _iris_codes;
	}
	if ( _fragility_maps )
	{
		for ( unsigned int i = 0; i < _nb_classes; ++ i )
		{
			cvReleaseImage ( _fragility_maps + i );
		}
		delete[] _fragility_maps;
	}
	if ( _fragility_rate )
	{
		for ( unsigned int i = 0; i < _nb_classes; ++ i )
		{
			if ( _fragility_rate[i] )
				delete[] _fragility_rate[i];

		}
		delete[] _fragility_rate;
	}
	if ( _iris_codes_bis)
	{
		for ( unsigned int i = 0; i < _nb_classes; ++ i )
		{
			if ( _iris_codes_bis[i]  )
				delete[] _iris_codes_bis[i];
		}
		delete[] _iris_codes_bis;
	}
			
	if ( _masks)
	{
		for ( unsigned int i = 0; i < _nb_classes; ++ i )
		{
			if ( _masks[i]  )
				delete[] _masks[i];
		}
		delete[] _masks;
		
	}
	
	if ( _fragile_bits )
	{
		for ( unsigned int i = 0; i < _nb_classes; ++ i )
		{
			if ( _fragile_bits[i]  )
				delete[] _fragile_bits[i];
		}
		delete[] _fragile_bits;
		
	}
	
	if ( _fragile_bit_thresholds )
		delete[] _fragile_bit_thresholds;
	
}

void c_database :: initialize()
{
	_nb_classes = 0;
	_class_names = NULL;
	_names = NULL;
	_fragile_bit_thresholds = 0;
	_previous_fragile_bit_rate = -1;
	_iris_codes = NULL;
	_fragility_maps = NULL;
	_fragility_rate = NULL;
	err_stream = NULL;
	_iris_codes_bis = NULL;
	_masks = NULL;
	_fragile_bits = NULL;
	
}

unsigned int c_database :: get_nb_classes(	const char * rep ) const
{
	stringstream stream;
	DIR * dir;
	struct dirent * ent;
	dir = opendir( rep );
	if ( ! dir )
	{
		cout << "Error: Unable to load " << rep << "!" << endl;
		return 1;
	}
	
	//Nombre
	unsigned int nb = 0;
	while ( ( ent = readdir( dir ) ) != NULL )
	{
		if ( strcmp( ent->d_name, ".." ) && strcmp( ent->d_name, "." ) )
		{
			nb ++;
		}
	}
	closedir(dir);
	return nb;
}

void c_database :: compute_fragility_rate( void )
{
	_fragility_rate[_nb_classes] = new double[256];
	hist.compute(	&( GET_IMAGE_PIXEL( _fragility_maps[_nb_classes], 0, 0, unsigned char ) ),
					GET_IMAGE_WIDTH( _iris_codes[_nb_classes] ),
					GET_IMAGE_HEIGHT( _fragility_maps[_nb_classes] ),
					GET_IMAGE_WIDTH_STEP( _fragility_maps[_nb_classes] ) );
	_fragility_rate[_nb_classes][0] = 0;
	_fragility_rate[_nb_classes][255] = 1;
	_fragility_rate[_nb_classes][1] = hist.data()[1];
	for ( unsigned int i = 2; i < 256; ++ i )
		_fragility_rate[_nb_classes][i] = hist.data()[i] + _fragility_rate[_nb_classes][i - 1];
	
	for ( unsigned int i = 1; i < 256; ++ i )
	{
		_fragility_rate[_nb_classes][i] /= _fragility_rate[_nb_classes][255];
		

	}



}



int c_database :: load (	const char * rep,
							const char * iris_code_file,
							const char * fragility_map_file )
{
	DIR * dir;
	struct dirent * ent;
	
	//Nombre de classes
	_nb_classes = get_nb_classes( rep );

	//Check ( Dossier vide )
	if ( _nb_classes == 0 )
		return 1;
	
	//Alloc mémoire
	alloc();

	//Nouveau parcours de la base
	dir = opendir( rep );

	_nb_classes = 0;
	while ( ( ent = readdir( dir ) ) != NULL )
	{
		
		if ( strcmp( ent->d_name, ".." ) && strcmp( ent->d_name, "." ) )
		{
			stringstream oss;
			oss << rep << "/" << ent->d_name ;
			//Test si c'est un répertoire
			if ( is_dir( oss.str().c_str() ) )
			{
				oss.str("");
				//Récupération du nom et du nom de la classe
				_names[_nb_classes] = ent->d_name;
				char buffer[8];
					memcpy ( buffer, ent->d_name, 7 );
					buffer[7] = '\0';
				_class_names[_nb_classes] = buffer;

				//Chargement de l'iris code
				oss << rep << "/" << ent->d_name << "/" << iris_code_file;
				_iris_codes[_nb_classes] = cvLoadImage( 	oss.str().c_str(),
															0 );
				oss.str("");
				
				//Chargement de la carte de fragilité
				oss << rep << "/" << ent->d_name << "/" << fragility_map_file;
				_fragility_maps[_nb_classes] = cvLoadImage( 	oss.str().c_str(),
																0 );
				oss.str("");

				//Check de la validité des deux
				if ( 	! _iris_codes[_nb_classes]  	||
						! _fragility_maps[_nb_classes] 	)
				{
					if ( _iris_codes[_nb_classes] )
					{
						cvReleaseImage( &_iris_codes[_nb_classes]  );
						_iris_codes[_nb_classes] = 0;
					}
					
					
					if ( _fragility_maps[_nb_classes] )
					{
						cvReleaseImage( &_fragility_maps[_nb_classes]  );
						_fragility_maps[_nb_classes] = 0;
					}	
				}
				else
				{
					compute_fragility_rate();

					
					
				}
				_nb_classes ++;
			}
		}
					
	}
	closedir(dir);
	return 0;
	
	
}

void c_database :: alloc()
{
	_class_names = new string[_nb_classes];
	_names = new string[_nb_classes];
	_iris_codes = new IplImage*[_nb_classes];
	_fragility_maps = new IplImage*[_nb_classes];	
	_fragility_rate = new double*[_nb_classes];
	_iris_codes_bis = new unsigned long long*[_nb_classes];
	_masks = new unsigned long long*[_nb_classes];
	_fragile_bits = new unsigned long long*[_nb_classes];
	
	
	
	for ( unsigned int i = 0; i < _nb_classes; ++ i )
	{
		_iris_codes_bis[i] = NULL;
		_masks[i] = NULL;
		_fragile_bits[i] = NULL;
		_fragility_rate[i] = new double[256];
	}
	
	
	_fragile_bit_thresholds = new double[_nb_classes];
}

double c_database :: get_fbt( double fbr, unsigned id )
{
	
	if ( _fragility_rate[id] == NULL )
		return 255;
	else
	{
		for ( unsigned int i = 1; i < 256; ++ i )
		{
			if ( _fragility_rate[id][i] >= fbr )
			{
				//Calcul des seuils
				if ( _iris_codes_bis[id] || _masks[id] || _fragile_bits[id] )
				{
					delete[] _iris_codes_bis[id];
					delete[] _masks[id];
					delete[] _fragile_bits[id];	
					_iris_codes_bis[id] = NULL;
					_masks[id] = NULL;
					_fragile_bits[id] = NULL;
					
					
				}
				
				if ( _iris_codes[id] && _fragility_maps[id] )
				{
					_iris_codes_bis[id] = new unsigned long long[ ((_iris_codes[id]->width + 63) /64 ) *  _iris_codes[id]->height];
					_masks[id] = new unsigned long long[ ((_iris_codes[id]->width + 63) /64 ) *  _iris_codes[id]->height];
					_fragile_bits[id] = new unsigned long long[ ((_iris_codes[id]->width + 63) /64 ) *  _iris_codes[id]->height];
					int q = distance::convert(	_iris_codes_bis[id],
												_masks[id],
												_fragile_bits[id],
												width_step,
												(const IplImage*) _iris_codes[id],
												(const IplImage*) _fragility_maps[id],
												i );
				}
				return i;
			}
		}
		
	}	
	return 255;
}

void c_database :: compute_fragile_bit_threshold( double fragile_bit_rate )
{

	if ( _previous_fragile_bit_rate != fragile_bit_rate )
	{
		for ( unsigned int i = 0; i < _nb_classes; ++ i )
		{

			_fragile_bit_thresholds[i] = get_fbt( fragile_bit_rate, i );
		}

		_previous_fragile_bit_rate = fragile_bit_rate;
	}
}

double c_database :: get_fragile_bit_threshold ( unsigned int id ) const
{
	if ( id >= _nb_classes )
		return -1;
	else
		return _fragile_bit_thresholds[id];
}
