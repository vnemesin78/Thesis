#include "lib_iris.hpp"
#include <dirent.h>
#include <ctime>
int main( int argc, char ** argv )
{
	c_get_iris_template_pthread obj;
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
		cout << "argv[1] : Metadata directory " << endl;
		cout << "argv[2] : Save directory " << endl;
		cout << "Author: Valérian Némesin." << endl;
		unsigned int id = 3;
		id = obj.help( cout, id );
		return 0; 
	}	
	
	if ( argc < 3 )
	{
		cout << "Error : missing argument(s)" << endl;
		return 1;
	}
	

	
	if ( obj.setup( argc - 3, argv + 3, argv[2], &cout, true ) )
		return 1;
		

	//Création du rép. de sauvegarde
	mkdir ( argv[2],
			014777 );
	
	//Chargement des paramères
	api_parameters params;
	for ( int i = 3; i < argc; ++ i )
		params.load(argv[i]);
	
	
	//Lecture du dossier et segmentation
	DIR * dir;
	struct dirent * ent;
	dir = opendir( argv[1] );
	
	if ( ! dir )
	{
		cout << "Error: Unable to load " << argv[1]  << "!" << endl;
		return 1;
	}
	
	unsigned int nb_videos = 0, n_video = 0;
	while ( ( ent = readdir( dir ) ) != NULL )
	{
		stringstream oss;
		oss << argv[1] << "/" << ent->d_name;
		if ( strcmp( ent->d_name, ".." ) && strcmp( ent->d_name, "." ) )
			nb_videos ++;
	}
	closedir(dir);
	dir = opendir( argv[1] );
	while ( ( ent = readdir( dir ) ) != NULL )
	{

		//Fichier
		stringstream oss;
		oss << argv[1] << "/" << ent->d_name;
		if ( strcmp( ent->d_name, ".." ) && strcmp( ent->d_name, "." ) )
		{
			n_video ++;
			system("clear");
			cout << "( " << n_video << " / " << nb_videos << " )" << endl; 
			
			api_parameters params;
			if ( ! params.load( oss.str().c_str() ))
			{
				string path, type;
				//Recup. du chemin et du type 
				if (	! api_get_string( params, "type", &type, &cout ) &&
						! api_get_string( params, "path", &path, &cout )	)
				{
					obj.test();
					//Récupération du dossier
					char * buffer = NULL;
					//Détection du dernier / dans path
					unsigned int 	pos_b = 0,
									pos_p = path.size();
					for ( unsigned int i = 0; i < path.size(); ++ i )
					{
						if ( path[i] == '/' )
							pos_b = i + 1;
						else if(  path[i] == '.' )
							pos_p = i;
					}
					
					
					
					if ( pos_b < pos_p )
					{
						buffer = new char[pos_p - pos_b + 1];
						memcpy( buffer, path.c_str() + pos_b, pos_p - pos_b  );
						//cerr << pos_p - pos_b << endl;
						buffer[pos_p - pos_b] = '\0';
						stringstream oss2;
						
						oss2 << argv[2] << "/" << type << '_' << buffer;
						mkdir ( oss2.str().c_str(),
									014777 );
						oss2.str("");		
						oss2 << argv[2] << "/" << type << '_' << buffer << "/log.txt";			
						ofstream log_file(oss2.str().c_str());
						obj.set_error_stream(log_file);		
						oss2.str("");
						
						
						
						
						//Segmentation
						while ( obj.run( path.c_str(), end ) );
						oss2 << argv[2] << "/" << type << '_' << buffer;
						obj.save( oss2.str().c_str() );
						oss2.str("");	
						delete[] buffer;
						log_file.close();
					}
					
					
				}
			}
			
			
		}
	}
	closedir(dir);
	return 0;
}
