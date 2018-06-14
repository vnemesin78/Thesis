#include "c_learning.hpp"
#include <cstring>
#include <ctime>
#include <dirent.h>
#include <sstream>
iris_class_data :: iris_class_data( 	
	const string & name,
	unsigned int nb_img_max )
{
	
	iris_class_data :: name = name;
	type = name[0];
	nb_images = 0;
	nb_images_max = nb_img_max;
	paths = new string[nb_img_max];
}

int iris_class_data :: add_image( 	const string & img_path )
{
	if ( nb_images == nb_images_max )
			return 1;
	
	paths[nb_images] = img_path;
	nb_images ++;
	
	return 0;
}

iris_class_data :: ~iris_class_data()
{
	if ( paths )
		delete[] paths;
	paths = 0;
	
}


c_learning :: c_learning ( 	)
{
	initialize();
}

c_learning :: c_learning ( 	
	int argc,
	char ** argv,
	unsigned int nb_classes_max,
	unsigned int nb_img_max )
{
	
	initialize();
	setup ( argc, 
			argv, 
			nb_classes_max, 
			nb_img_max );
}

unsigned int c_learning :: help( 	ostream & out, 
										unsigned int id ) const
{
	out << "Argv[" << id ++ << "] : metadata folder" << endl;
	out << "Argv[" << id ++ << "] : save folder" << endl;
	
	id = obj.help( out, id );
	
	return id;
}

int c_learning :: setup( 	
	int argc,
	char ** argv,
	unsigned int nb_classes_max,
	unsigned int nb_img_max )
{
	free();
	initialize();
	
	_nb_img_max = nb_img_max;
	_nb_classes_max = nb_classes_max;
	data = new iris_class_data*[nb_classes_max];
	
	if ( argc < 3 )
	{
		cout << "Error: Missing argument(s)!" << endl;
		return 1;
	}
	
	//Création du rép. de sauvegarde
	mkdir (argv[2], 014777 );
	_save_rep = argv[2];
	//Fichier de logs
	{
		stringstream oss;
		oss << argv[2] << "/log.txt";
		log = new ofstream( oss.str().c_str() );
	}

	//Setup de l'obj de seg.
	if (	obj.setup( 	argv + 2, 
						argc - 2, 
						log ) ) 
	{
		stringstream oss;
		oss << argv[2] << "/log.txt";
		cout << "Error: Check " << oss.str().c_str() << " for more information!" << endl;
		return 1;
	}
	
	
	//Chargemnent des métadata
	{
		DIR * dir = opendir( argv[1] );
		if ( ! dir )
		{
			cout << "Error: Can't open metadata dir "<< argv[1] << " !" << endl;
			return 1;
			
		} 
		struct dirent * ent;
		while ( ( ent = readdir( dir ) ) != NULL )
		{
			if ( strcmp( ent->d_name, ".." ) && strcmp( ent->d_name, "." ) )
			{
				api_parameters params;
				stringstream oss;
				oss << argv[1] << "/" << ent->d_name;

				if ( ! params.load(oss.str().c_str() ) )
				{

					string c_name, 
							path;
					if (	! api_get_string( params, "path", &path, &cout )		&&
							! api_get_string( params, "class", &c_name, &cout )	)
					{
						
						bool add = false;
						for ( unsigned int i = 0; i < _nb_classes; ++ i )
						{
							if ( data[i]->name == c_name )
							{
								add = true;
								if ( data[i]->add_image( path ) )
								{
									cout << c_name << endl;
									cout << "Error! Check iris_class_data parameters!" << endl;
								}
								break;
								
							}
						}
						if ( ! add )
						{
							if ( _nb_classes == _nb_classes_max )
							{
									cout << "Error! Nb max classes "<< _nb_classes_max<<" !" << endl;
								
							}
							else
							{
								data[_nb_classes] = 
									new iris_class_data( 	c_name, 
															_nb_img_max );
								data[_nb_classes]->add_image( path );
								_nb_classes ++;

								
								
							}
						}
						
						
					}
				}
			}
		}
		
		
		closedir (dir);
	}
	
	
	return 0;
}

int c_learning :: segment( )
{
	unsigned int total = 0, nb_fails = 0;
	stringstream oss;
	oss << _save_rep << "/" << "seg.txt";
	ofstream seg_data(oss.str().c_str());
	oss.str("");
	
	oss << _save_rep << "/" << "time.txt";
	ofstream time_file(oss.str().c_str());
	oss.str("");
	
	
	
	//Segmentation par classes
	for ( unsigned int i = 0; i < _nb_classes; ++ i )
	{
		unsigned int success = 0;
		oss << _save_rep << "/" << data[i]->name;
			mkdir( oss.str().c_str(), 014777 );
		oss.str("");
		oss << _save_rep << "/" << data[i]->name << "/" << "log.txt";
			ofstream log_class( oss.str().c_str() );
		
		oss.str("");
		


		seg_data << "#" << data[i]->name << endl;
		seg_data << data[i]->name << "::nb_images = " << data[i]->nb_images << endl;
	

		obj.set_error_stream( log_class );

		oss << _save_rep << "/" << data[i]->name;

		for ( unsigned int j = 0; j < data[i]->nb_images; ++ j )
		{
			
			
			
			
			log_class << ((*data[i])[j])->c_str() << "\t";
			clock_t time_1 = clock();
			IplImage * image = 
				cvLoadImage( ((*data[i])[j])->c_str(), 0 );
			if ( image ) 
			{

				log_class << "OPEN_OK" << "\t";
				if ( ! obj.segment( image, 
									((*data[i])[j])->c_str() ) )
				{
					log_class << "SEG_OK" << "\t";		
					success ++;
				}
				else
				{
					nb_fails ++;
					log_class << "SEG_FAILED" << "\t";	
				}
				obj.data().save( oss.str().c_str(), j );	
				
				clock_t time_2 = clock();	
				
				time_file << "times(" << total + j + 1 << ") = " << (time_2 - time_1) / ( (double) CLOCKS_PER_SEC ) << ";" << endl;
				
				cvReleaseImage( &image );
			}
			else
			{
				log_class << "OPEN_FAILED" << "\t";				
				
			}
			log_class << endl;
		
			system("clear");
			cout << " ( " << i + 1 << " / " << _nb_classes << " ) " << endl;
			total ++;
			if ( total > 0 )
			{
				double f = ( (10000 * nb_fails) / total ) / 100.0;
				cout << " fails: " << f << "%" << endl;
				if ( f < 2.5 )
				{
					cout <<
	"                    __ooooooooo__" << endl <<
	"                 oOOOOOOOOOOOOOOOOOOOOOo" << endl <<
	"             oOOOOOOOOOOOOOOOOOOOOOOOOOOOOOo" << endl <<
	"          oOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOo" << endl <<
	"        oOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOo" << endl <<
	"      oOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOo" << endl <<
	"     oOOOOOOOOOOO*  *OOOOOOOOOOOOOO*  *OOOOOOOOOOOOo" << endl <<
	"    oOOOOOOOOOOO      OOOOOOOOOOOO      OOOOOOOOOOOOo" << endl <<
	"    oOOOOOOOOOOOOo  oOOOOOOOOOOOOOOo  oOOOOOOOOOOOOOo" << endl <<
	"   oOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOo" << endl <<
	"   oOOOO     OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO     OOOOo" << endl <<
	"   oOOOOOO OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OOOOOOo" << endl <<
	"    *OOOOO  OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO  OOOOO*" << endl <<
	"    *OOOOOO  *OOOOOOOOOOOOOOOOOOOOOOOOOOOOO*  OOOOOO*" << endl <<
	"     *OOOOOO  *OOOOOOOOOOOOOOOOOOOOOOOOOOO*  OOOOOO*" << endl <<
	"      *OOOOOOo  *OOOOOOOOOOOOOOOOOOOOOOO*  oOOOOOO*" << endl <<
	"       *OOOOOOOo  *OOOOOOOOOOOOOOOOO*  oOOOOOOO*" << endl <<
	"         *OOOOOOOOo  *OOOOOOOOOOO*  oOOOOOOOO*   " << endl <<
	"             *OOOOOOOOo           oOOOOOOOO*   " << endl <<
	"                 *OOOOOOOOOOOOOOOOOOOOO* " << endl <<
	"                      \"\"ooooooooo\"" << endl;
				}
				else if ( f < 5 )
				{
					
					cout <<		
	"                    @@@" << endl <<
	"                             @@@" << endl <<
	"                              @@@" << endl <<
	"                              @@@" << endl <<
	"                      @@@@@@@@@@@@@@@@@@@@@@" << endl <<
	"                    @@@@@@@@@@@@@@@@@@@@@@@@@@" << endl <<
	"                  @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@" << endl <<
	"                @@@@@@@@ @@@@@@@@@@@@@@@@ @@@@@@@@" << endl <<
	"              @@@@@@@@@   @@@@@@@@@@@@@@   @@@@@@@@@" << endl <<
	"            @@@@@@@@@@     @@@@@@@@@@@@     @@@@@@@@@@" << endl <<
	"           @@@@@@@@@@       @@@@  @@@@       @@@@@@@@@@" << endl <<
	"           @@@@@@@@@@@@@@@@@@@@    @@@@@@@@@@@@@@@@@@@@" << endl <<
	"           @@@@@@@@@@@@@@@@@@        @@@@@@@@@@@@@@@@@@" << endl <<
	"           @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@" << endl <<
	"           @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@" << endl <<
	"           @@@@@@@@@ @@@@@@@@@@@@@@@@@@@@@@@@ @@@@@@@@@" << endl <<
	"            @@@@@@@@  @@ @@ @@ @@ @@ @@ @@ @  @@@@@@@@" << endl <<
	"              @@@@@@@                        @@@@@@@" << endl <<
	"                @@@@@@  @@ @@ @@ @@ @@ @@ @ @@@@@@" << endl <<
	"                  @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@" << endl <<
	"                    @@@@@@@@@@@@@@@@@@@@@@@@@@" << endl <<
	"                      @@@@@@@@@@@@@@@@@@@@@@" << endl;
					
				}
				else
				{
	cout << 
	"             ____"<< endl <<
	"      _,-ddd888888bbb-._"<< endl <<
	"    d88888888888888888888b"<< endl <<
	"  d888888888888888888888888b"<< endl <<
	" 6888888888888888888888888889"<< endl <<
	" 68888b8""8q8888888p8""8d88889"<< endl <<
	" `d8887     p88888q     4888b\'"<< endl <<
	"  `d8887    p88888q    4888b\'"<< endl <<
	"    `d887   p88888q   488b\'"<< endl <<
	"      `d8bod8888888dob8b\'"<< endl <<
	"        `d88888888888d\'"<< endl <<
	"          `d8888888b\' hjw"<< endl <<
	"            `d8888b\' `97"<< endl <<
	"              `bd\'"<< endl;
				}
		
			}
		}
		log_class.close();
		oss.str("");
		seg_data << data[i]->name << "::nb_success = " << success << endl;
	}
	seg_data << "#Global" << endl;
	seg_data << "nb_images = " << total << endl;
	seg_data << "nb_fails = " << nb_fails << endl;
	seg_data << "nb_success = " << total - nb_fails << endl;
	time_file.close();
	seg_data.close();
	return 0;
}

c_learning :: ~c_learning()
{
	free();
	initialize();
}

void c_learning :: free()
{
	if (data)
	{
		for (unsigned int i = 0; i < _nb_classes; ++ i )
		{
			delete data[i];
		}
		delete[] data;
	}
	
	if ( log )
	{
		log->close();
		delete log;
	}		
}

void c_learning :: initialize()
{
	data = 0;
	log = 0;
	_nb_classes = 0;
	_nb_classes_max = 0;
	_nb_img_max = 0;
	
}


