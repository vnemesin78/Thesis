#include "c_recognition.hpp"
#include <ctime>
c_recognition :: c_recognition()
{
	initialize();
}

c_recognition :: c_recognition( int argc, 
								char ** argv,
								const char * iris_code_fn,
								const char * fragility_map_fn,
								unsigned int nb_fused_images,
								unsigned int nb_thresholds  )
{
	initialize();
	setup( 	argc,
			argv,
			iris_code_fn,
			fragility_map_fn,
			nb_fused_images,	
			nb_thresholds );
}

int c_recognition :: setup(  	int argc, 
								char ** argv,
								const char * iris_code_fn,
								const char * fragility_map_fn,
								unsigned int nb_fused_images,
								unsigned int nb_thresholds  )
{

	free();

	initialize();
	
	videos = new c_database;
	images = new c_matching;
	
	//Chargement des paramètres
	api_parameters params;
	for ( int i = 1; i < argc; ++ i )
		params.load(argv[i]);

	//videos
	{
		stringstream oss, oss2;
		if ( nb_fused_images != 0 )
		{
			oss << iris_code_fn << "_" << nb_fused_images << ".png";
			oss2 << fragility_map_fn << "_" << nb_fused_images << ".png";
		}
		else
		{
			oss << iris_code_fn << ".png";
			oss2 << fragility_map_fn << ".png";	
		}
		
		if ( videos->setup( argv[2], 
							oss.str().c_str(), 
							oss2.str().c_str() ) )
		{

			delete videos;
			delete images;
			return 1;
		}
		

		if ( images->setup( argv[1], 
							oss.str().c_str(), 
							oss2.str().c_str(),
							params,
							&cout   ) )
		{

			delete videos;
			delete images;
			return 1;
		}

	}
	
	_id = nb_fused_images;
	_nb_thresholds = nb_thresholds;
	
	truth = new unsigned int[videos->nb_classes()];
	orders = new unsigned int[videos->nb_classes() * images->nb_classes()];
	top_rank = new unsigned int[videos->nb_classes() * images->nb_classes()];
	distances = new double[videos->nb_classes() * images->nb_classes()];
	
	v_p = new double[nb_thresholds];
	v_n = new double[nb_thresholds];
	f_p = new double[nb_thresholds];
	f_n = new double[nb_thresholds];
	
	thresholds = new double[nb_thresholds];
	i_distances = new double[ videos->nb_classes() * (images->nb_classes() - 1)];
	g_distances = new double[ videos->nb_classes()];
	
	//Construction de la vérité terrain
	for ( unsigned int i = 0; i < videos->nb_classes(); ++ i )	
	{	
		unsigned int j;
		for ( j = 0; j < images->nb_classes(); ++ j )
		{
			if ( !strcmp( videos->class_name(i), images->class_name(j) ) )
				break;		
		}

		truth[i] = j; 
	}

	return 0;
}

int c_recognition :: match( double fbr, 
							double alpha)
{
	double matchings = 0;
	double time = 0;
	//Calcul des FBR
	if ( use_fbr() )
	{
		videos->compute_fragile_bit_threshold( fbr );
		images->compute_fragile_bit_threshold( fbr );
	}

	//Alpha
	if ( use_alpha() )
		images->set_alpha( alpha );
	
	
	//Matching
	for ( _nb_video = 0;  _nb_video < videos->nb_classes(); ++ _nb_video )	
	{
		//~ system("clear");
		clock_t time_1 = clock();
		unsigned int q_truth;
		if ( images->get_distance() == (distance::function_prototype_bis) & (distance :: Hamming_expectation) )
		{
			images->matching( 	distances + _nb_video * images->nb_classes(),
								orders + _nb_video * images->nb_classes(),
								q_truth,
								videos->iris_code(_nb_video), 
								videos->fragility_map(_nb_video), 
								videos->get_fragile_bit_threshold(_nb_video),
								videos->class_name(_nb_video),
								videos->name(_nb_video)
								);
		}
		else
		{		
			images->matching_opt( 	distances + _nb_video * images->nb_classes(),
									orders + _nb_video * images->nb_classes(),
									q_truth,
									videos->iris_code_bis(_nb_video), 
									videos->mask(_nb_video), 
									videos->fragile_bits(_nb_video),
									videos->class_name(_nb_video),
									videos->name(_nb_video)
									);				
								
		}
		clock_t time_2 = clock();
		matchings += images->nb_classes();
		time += (time_2 - time_1) / ( (double) CLOCKS_PER_SEC );
		
		
	}
	//Roc
	sort_distances();

	get_roc_data();

	//Top rank
	get_top_rank();


	return 0;
}

		
int c_recognition :: save( const char * rep ) const
{
	unsigned int id = _id;
	stringstream poulet_terminator;
	FILE * file = NULL;
	
	poulet_terminator << rep << "/" << "data_nb_img_" << id << ".m";
	
	file = fopen( poulet_terminator.str().c_str(), "w" );
	if ( !file )
		return 1;
	poulet_terminator.str("");
	
	
	
	//Sauvegarde du tableau des distances
	{
		gsl_matrix toto;
		toto.size1 = videos->nb_classes();
		toto.size2 = images->nb_classes();
		toto.data = distances;
		toto.tda = images->nb_classes();
		
		poulet_terminator << "distances_" << id;
		fprintf( file, "%s = ", poulet_terminator.str().c_str() );
		gsl_matrix_fprintf_(	file,
								"[^];\n^\t^\n^%lf",
								&toto );
		poulet_terminator.str("");
	}

	//Sauvegarde des ordres
	{
		gsl_matrix * ord_mat = gsl_matrix_alloc( 	videos->nb_classes(), 
													images->nb_classes());
		
		for ( unsigned int i = 0; i < videos->nb_classes(); ++ i )
		{
			for ( unsigned int j = 0; j < images->nb_classes(); ++ j )
			{
				ord_mat->data[ ord_mat->tda * i + j] = orders[ images->nb_classes() * i + j] + 1;
			}
			
		}
		
		poulet_terminator << "orders_" << id;
		fprintf( file, "%s = ", poulet_terminator.str().c_str() );
		gsl_matrix_fprintf_(	file,
								"[^];\n^\t^\n^%lf",
								ord_mat );
		poulet_terminator.str("");
		
		
		gsl_matrix_free(ord_mat);
	}

	//Sauvegarde de la vérité terrain
	{
		gsl_matrix * ord_mat = gsl_matrix_alloc( 	1, 
													videos->nb_classes());
		
		for ( unsigned int i = 0; i < videos->nb_classes(); ++ i )
		{
			ord_mat->data[i] = truth[i] + 1;
		}
		
		fprintf( file, "%s = ", "truth" );
		gsl_matrix_fprintf_(	file,
								"[^];\n^\t^\n^%lf",
								ord_mat );
		gsl_matrix_free(ord_mat);
	}

	//Sauvegarde du TOP RANK
	{
		gsl_matrix * ord_mat = gsl_matrix_alloc( 	1, 
													videos->nb_classes());
		
		for ( unsigned int i = 0; i < videos->nb_classes(); ++ i )
		{
			ord_mat->data[i] = top_rank[i] + 1;
		}
		
		poulet_terminator << "top_rank_" << id;
		fprintf( file, "%s = ", poulet_terminator.str().c_str() );
		gsl_matrix_fprintf_(	file,
								"[^];\n^\t^\n^%lf",
								ord_mat );
		poulet_terminator.str("");		
		gsl_matrix_free(ord_mat);
	}

	//Sauvegarde des seuils
	{
		gsl_matrix toto;
		toto.size1 = 1;
		toto.size2 = _nb_thresholds;
		toto.data = thresholds;
		toto.tda = _nb_thresholds;
		
		poulet_terminator << "s";
		fprintf( file, "%s = ", poulet_terminator.str().c_str() );
		gsl_matrix_fprintf_(	file,
								"[^];\n^\t^\n^%lf",
								&toto );
		poulet_terminator.str("");
	}
	//Sauvegarde des seuils
	{
		gsl_matrix toto;
		toto.size1 = 1;
		toto.size2 = _nb_thresholds;
		toto.data = v_p;
		toto.tda = _nb_thresholds;
		
		poulet_terminator << "v_p_" << id;
		fprintf( file, "%s = ", poulet_terminator.str().c_str() );
		gsl_matrix_fprintf_(	file,
								"[^];\n^\t^\n^%lf",
								&toto );
		poulet_terminator.str("");
	}
	//Sauvegarde des seuils
	{
		gsl_matrix toto;
		toto.size1 = 1;
		toto.size2 = _nb_thresholds;
		toto.data = v_n;
		toto.tda = _nb_thresholds;
		
		poulet_terminator << "v_n_" << id;
		fprintf( file, "%s = ", poulet_terminator.str().c_str() );
		gsl_matrix_fprintf_(	file,
								"[^];\n^\t^\n^%lf",
								&toto );
		poulet_terminator.str("");
	}
	//Sauvegarde des seuils
	{
		gsl_matrix toto;
		toto.size1 = 1;
		toto.size2 = _nb_thresholds;
		toto.data = f_n;
		toto.tda = _nb_thresholds;
		
		poulet_terminator << "f_n_" << id;
		fprintf( file, "%s = ", poulet_terminator.str().c_str() );
		gsl_matrix_fprintf_(	file,
								"[^];\n^\t^\n^%lf",
								&toto );
		poulet_terminator.str("");
	}
	
	//Sauvegarde des seuils
	{
		gsl_matrix toto;
		toto.size1 = 1;
		toto.size2 = _nb_thresholds;
		toto.data = f_p;
		toto.tda = _nb_thresholds;
		
		poulet_terminator << "f_p_" << id;
		fprintf( file, "%s = ", poulet_terminator.str().c_str() );
		gsl_matrix_fprintf_(	file,
								"[^];\n^\t^\n^%lf",
								&toto );
		poulet_terminator.str("");
	}

	//Sauvegarde des seuils
	{
		gsl_matrix toto;
		toto.size1 = 1;
		toto.size2 = videos->nb_classes() * ( images->nb_classes() - 1 );
		toto.data = i_distances;
		toto.tda = videos->nb_classes() * ( images->nb_classes() - 1 );
		
		poulet_terminator << "i_distances_" << id;
		fprintf( file, "%s = ", poulet_terminator.str().c_str() );
		gsl_matrix_fprintf_(	file,
								"[^];\n^\t^\n^%lf",
								&toto );
		poulet_terminator.str("");
	}
	
	{
		gsl_matrix toto;
		toto.size1 = 1;
		toto.size2 = videos->nb_classes();
		toto.data = g_distances;
		toto.tda = videos->nb_classes();
		
		poulet_terminator << "g_distances_" << id;
		fprintf( file, "%s = ", poulet_terminator.str().c_str() );
		gsl_matrix_fprintf_(	file,
								"[^];\n^\t^\n^%lf",
								&toto );
		poulet_terminator.str("");
	}
	fprintf( file, "err = %lf;\n", err );
	
	//Sauvegarde des noms
	{
		fprintf( file, "video_names = {};" );
		for ( unsigned int i = 0; i < videos->nb_classes(); ++ i )
		{
			fprintf( 	file, 
						"video_names(%u) = \"%s\";\n", 
						i + 1, 
						videos->name(i) );
		
		}
	}
	
	//Sauvegarde des noms
	{
		fprintf( file, "class_names = {};" );
		for ( unsigned int i = 0; i < images->nb_classes(); ++ i )
		{
			fprintf( 	file, 
						"class_names(%u) = \"%s\";\n", 
						i + 1, 
						images->name(i) );
		}
	}
	
	
	fclose(file);
	
	return 0;
}

c_recognition :: ~c_recognition()
{
	free();
	initialize();

}

void c_recognition :: initialize()
{
	videos = 0;
	images = 0;
	_nb_video = 0;		
	_id = 0;
	_nb_thresholds = 0;
	truth = 0;
	orders = 0;
	top_rank = 0;
							
	distances = 0;
	v_p = 0;
	f_p = 0;
	v_n = 0;
	f_n = 0;
	thresholds = 0;
	i_distances = 0;
	g_distances = 0;
	
}

void c_recognition :: free()
{
	if ( videos )
		delete videos;
	
	if ( images )
		delete images;
		
	if ( truth )
		delete[] truth;
		
	if ( orders )
		delete[] orders;
		
	if ( top_rank )
		delete[] top_rank;
		
	if ( distances )
		delete[] distances;
	
	if ( v_p )
		delete[] v_p;
	
	if ( v_n )
		delete[] v_n;
		
	if ( f_p )
		delete[] f_p;
		
	if ( f_n )
		delete[] f_n;
		
	if ( i_distances )
		delete[] i_distances;
		
	if ( g_distances )
		delete[] g_distances;
	
	if ( thresholds )
		delete[] thresholds;

}

void c_recognition :: sort_distances()
{
	unsigned int i1 = 0, i2 = 0;
	for (unsigned int i = 0; i < videos->nb_classes(); ++ i )
	{
		for ( unsigned int j = 0; j < images->nb_classes(); ++ j )
		{
			if ( truth[i] == j )
			{
				g_distances[i1] = distances[i * images->nb_classes() + j];				
				i1 ++;
			}
			else
			{
				i_distances[i2] = distances[i * images->nb_classes() + j];	
				i2 ++;
			}
			
		}
		
	}
}

void c_recognition :: get_roc_data()
{
	unsigned int nb_videos = videos->nb_classes();
	//Erreurs d'enrollement
	for ( unsigned int i = 0; i < videos->nb_classes(); ++ i ) 
	{
		if ( g_distances[i] == 1 )
		{
			nb_videos --;
		}
		
	}
	err = 1.0 - nb_videos / ( (double) videos->nb_classes() );
	unsigned int t_n = (images->nb_classes() - 1);

	
	
	for ( unsigned int i = 0; i <_nb_thresholds; ++ i )
	{
		unsigned int nb_i = 0, nb_g = 0;
		
		thresholds[i] = i / ( (double) _nb_thresholds );  
		
		for ( unsigned int j = 0; j < videos->nb_classes(); ++j )
		{
			if ( g_distances[j] != 1 )
			{
				if ( g_distances[j] < thresholds[i] )
					++ nb_g;
			}
			for ( unsigned int k = 1; k < t_n; ++k )
			{
				if ( i_distances[j * t_n + k] < thresholds[i] )
					++ nb_i;
			}
		}

		v_p[i] = ( (double) nb_g) / nb_videos;
		f_p[i] = ( (double) nb_i) / (t_n * nb_videos);
		v_n[i] = 1.0 - v_p[i];
		f_n[i] = 1.0 - f_p[i];
	}
}

void c_recognition :: get_top_rank()
{
	for (unsigned int i = 0; i < videos->nb_classes(); ++ i )
	{
		top_rank[i] = images->nb_classes();	
		for ( unsigned int j = 0; j < images->nb_classes(); ++ j )
		{
			if ( truth[i] == orders[i * images->nb_classes() + j] )
			{
				top_rank[i] = j;				
				break;
			}
		}
	}
}

unsigned int c_recognition ::  help ( ostream & out, unsigned int id ) const
{
	out << "argv[" << id ++ << "- N] : .cfg files " << endl;
	return id;
	
}
