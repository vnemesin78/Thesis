#include "lib_iris.hpp"
#include "lib_api.hpp"
#include "lib_image.hpp"
#include <cstring>

#ifdef __linux__
	#define PATH_SEP '/'
#else
	#define PATH_SEP '\\'
#endif

int get_directories( char *** dir_names,
					 unsigned int & nb_dir,
					 const char * path, 
					 char sep = PATH_SEP )
{
	unsigned int nb = 0,
				 size = strlen( path );
	
	for ( unsigned int i = 0; i < size; ++ i )
		if ( path[i] == sep )
			if ( path[i + 1] != sep )
				nb ++;
	
	//Alloc mémoire
	(*dir_names) = new char*[nb + 1];
	
	//Recherche des noms
	{
		unsigned int pos = 0;
		for ( unsigned int i = 0; i < nb; ++ i )
		{
			unsigned int p_pos = pos;
			pos = gsl_matrix_io_split( 	path + pos, 
										sep);
			pos += p_pos;
			(*dir_names)[i] = new char[ pos - p_pos + 1];
			memcpy( (*dir_names)[i], 
					path + p_pos, 
					pos - p_pos + 1 );
			(*dir_names)[i][pos - p_pos] = '\0';
			pos ++;
		}
		(*dir_names)[nb] = new char[ size - pos + 1];
		memcpy( (*dir_names)[nb], path + pos, size - pos );
		(*dir_names)[nb][ size - pos ] = '\0';
	}
	
	//Suppression des /
	for ( unsigned int i = 0; i <= nb; ++ i )
	{
		for ( unsigned int j = 0; (*dir_names)[i][j] != '\0'; ++ j )
		{
			if ( (*dir_names)[i][j] == sep )
				(*dir_names)[i][j] = '\0';
		}
	}
	
	//Suppression des reps. vides
#ifdef __linux__
	if ( strlen( (*dir_names)[0] ) == 0 )
	{
		delete[] (*dir_names)[0];
		(*dir_names)[0] = new char[2];
		(*dir_names)[0][0] = sep;
		(*dir_names)[0][1] = '\0';
	}
#endif
	if ( strlen( (*dir_names)[nb] ) == 0  )
	{
		delete[] (*dir_names)[nb];
		nb --;
	}
	nb_dir = nb + 1;
	return 0;
}

int main ( int argc, char ** argv )
{
	char ** dir_names = NULL;
	unsigned int nb_dir = 0;
	
	if ( argc > 1 )
		if ( ! strcmp( argv[1], "-H" ) || ! strcmp( argv[1], "--help" ) )
		{
cout <<
"                                 HELP                                 \n"
" Arg[1] : Iris filename                                               \n"
" Arg[2] : Iris template filename                                      \n"
" Arg[3] : Seg. status filename                                        \n"
"\n"
" Brief\n"
" Enroll iris\n"
"Author: Valérian Némesin." << endl;
			return 0;
		}
		
	if ( argc < 4 )
	{
		cout << "Error: Missing arguments" << endl;
		return 1;
	}
	
	get_directories( &dir_names,
					 nb_dir,
					 argv[1] );
	int learning = 0;
	for ( unsigned int i = 0; i < nb_dir; ++ i )
	{
		if ( !strcmp(dir_names[i], "CASIA-Iris-Thousand" ) )
			learning = 1;
	}
	c_get_iris_template obj;
	if ( obj.default_setup(&cout) )
		return 1;
	//Matching
	if ( learning  == 0 )
	{

		//Chargement de l'image
		IplImage * image = cvLoadImage( argv[1], CV_LOAD_IMAGE_GRAYSCALE );
		if ( ! image )
		{
			cout << "Error : Can't open image " << argv[1] << endl;
			return 1;	
		}
		//Segmentation de l'image
		obj.segment( image, argv[1] );
		cvReleaseImage(&image);
		
		//Sauvegarde des résultats
		obj.safe_save( argv[2], argv[3], argv[1] );
		//mkdir("test",014777);
		//obj.data().save( "test", 0 );
	}
	
	//Apprentissage
	else if ( learning == 1 )
	{
		int n = 0;
		iris_data * i_data = new iris_data[10];
		unsigned int nb = 0;
		
		n += dir_names[nb_dir - 1][strlen(dir_names[nb_dir - 1]) - 6] - '0';
		n = n * 10 + dir_names[nb_dir - 1][strlen(dir_names[nb_dir - 1]) - 5] - '0';
		
		c_fusion_iris_template fusion;
		if ( fusion.default_setup() )
			return 1;

		stringstream oss;
		for ( unsigned int j = 0; j < nb_dir - 1; ++ j )
		{
			oss << dir_names[j] << "/";
		}
		DIR * dir;
		dir = opendir( oss.str().c_str() );
		if( !dir)
			return 1;
		struct dirent * ent;
		while ( ( ent = readdir( dir ) ) != NULL )
		{		
			if ( strcmp( ent->d_name, ".." ) && strcmp( ent->d_name, "." ) && strlen( ent->d_name) > 4 )
			{
				if ( !strcmp( ent->d_name + strlen(ent->d_name) - 4, ".jpg" ) )
				{
					stringstream oss2;
					oss2 << oss.str().c_str() << "/" << ent->d_name;
					
					//Chargement de l'image
					IplImage * image = cvLoadImage( oss2.str().c_str(), CV_LOAD_IMAGE_GRAYSCALE );
					if ( ! image )
					{
						cout << "Error : Can't open image " << argv[1] << endl;
						return 1;	
					}
					//Segmentation de l'image
					if ( ! obj.segment( image, argv[1] ) )
					{
						i_data[nb].setup( obj.data() );
						nb ++;
					}
					
					
					
					
					cvReleaseImage(&image);									
				}
			}
		}	
		closedir(dir);
		if ( fusion.fusion( i_data, nb ) )
		{
			ofstream f(argv[3]);
			f << 0;
			f.close();
		}
		else
		{
			//~ fusion.safe_save(argv[2], n, argv[1] );
			ofstream f(argv[3]);
			f << 1;
			f.close();
			//fusion.save("test", "fusion", "fusion_mask");
		}
		
		
		
		oss.str("");

	}
	
	for ( unsigned int i = 0; i < nb_dir; ++ i )
	{
		delete[] dir_names[i];
	}
	delete[] dir_names;
	
	
	return 0;
}
