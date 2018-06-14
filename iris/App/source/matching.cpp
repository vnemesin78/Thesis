#define NB_PROCESS 2
#include "../iris/lib_iris.hpp"

/**@struct
 * @brief
 * Gestion des paramètres du matching.
 * 
 **/
struct matching_parameters
{
	/**@fn
	 * @brief
	 * Constructeur ( Mise à zéro de tous les attributs).
	 * 
	 **/
	matching_parameters ( )
	{
		nb_thresholds = 0;
		alpha_min = 0;
		alpha_max = 0;
		nb_alpha = 0;
		fbr_min = 0;
		fbr_max = 0;
		nb_process = 1;
		nb_fbr = 0;
	} 
	
	
	/***@fn
	 * @brief
	 * Constructeur de l'objet.
	 * 
	 **/
	int setup ( api_parameters & params,
				const char * name_space = "matching",
				ostream * err_stream = &cout ) 
	{
		int q = 0;
		stringstream oss;
		oss << name_space << "::distance";
		if ( api_get_string ( params, 
							  oss.str().c_str(), 
							  &distance,
							  err_stream ) )
			q = 1;
		oss.str("");
		oss << name_space << "::fragility_map_filename";
		if ( api_get_string ( params, 
							  oss.str().c_str(), 
							  &fragility_map_fn,
							  err_stream ) )
			q = 1;
		oss.str("");
		oss << name_space << "::iris_code_filename";
		if ( api_get_string ( params, 
							  oss.str().c_str(), 
							  &iris_code_fn,
							  err_stream ) )
			q = 1;
		oss.str("");
		
		//NB_thresholds
		oss << name_space << "::nb_thresholds";
		if ( api_get_positive_integer (	params, 
										oss.str().c_str(), 
										&nb_thresholds,
										err_stream ) )
			q = 1;
		oss.str("");
		
		oss << name_space << "::nb_fused_images";
		if ( api_get_positive_integer (	params, 
										oss.str().c_str(), 
										&nb_fused_images) )
		{
			nb_fused_images = 0;
		}
		oss.str("");
		
		
		oss << name_space << "::nb_process";
		if ( api_get_positive_integer (	params, 
										oss.str().c_str(), 
										&nb_process,
										NULL ) )
		{
			if ( err_stream )
				*err_stream << "Warning : nb_process was set to 1." << endl;
			nb_process = 1;
		}
		oss.str("");
	
		
		
		
		//FBR
		if (	distance == "Hamming" 		|| 
				distance == "FBD" 			|| 
				distance == "Hamming_FBD"	)
		{
			oss << name_space << "::FBR_min";
			if ( api_get_double ( params, 
								  oss.str().c_str(), 
								  &fbr_min,
								  err_stream ) )
				q = 1;
			oss.str("");
			oss << name_space << "::FBR_max";
			if ( api_get_double ( params, 
								  oss.str().c_str(), 
								  &fbr_max,
								  err_stream ) )
				q = 1;
			oss.str("");
			oss << name_space << "::nb_FBR";
			if ( api_get_positive_integer (	params, 
											oss.str().c_str(), 
											&nb_fbr,
											err_stream ) )
				q = 1;
			oss.str("");
		}
		else
		{
			nb_fbr = 1;
			fbr_min = 0;
			fbr_max = 0;
		}
		
		//ALPHA
		if (	distance == "Hamming_FBD"	)
		{
			oss << name_space << "::alpha_min";
			if ( api_get_double ( params, 
								  oss.str().c_str(), 
								  &alpha_min,
								  err_stream ) )
				q = 1;
			oss.str("");
			oss << name_space << "::alpha_max";
			if ( api_get_double ( params, 
								  oss.str().c_str(), 
								  &alpha_max,
								  err_stream ) )
				q = 1;
			oss.str("");
			oss << name_space << "::nb_alpha";
			if ( api_get_positive_integer (	params, 
											oss.str().c_str(), 
											&nb_alpha,
											err_stream ) )
				q = 1;
			oss.str("");
			
		}
		else
		{
			nb_alpha = 1;
			alpha_min = 0;
			alpha_max = 0;
		}

		if ( q == 1 )
			return 1;
		return 0;
		
	}
	
	
	string iris_code_fn, 
		   fragility_map_fn,
		   distance;
		   
	unsigned int nb_fused_images;
	
	unsigned int nb_thresholds;		
	
	double alpha_min,
		   alpha_max;
	unsigned int nb_alpha;
	
	double fbr_min,
		   fbr_max;
	unsigned int nb_fbr;
	
	unsigned int nb_process;
};

struct matching_process
{
	/**@fn
	 * @brief
	 * Constructeur
	 * 
	 * 
	 */
	matching_process( )
	{
		obj = NULL;
	}
	/**@fn
	 * @param params : paramètres
	 * @param argc : nombre d'arguments
	 * @param argv : arguments
	 * @brief
	 * Setup
	 * 
	 * 
	 */
	int setup (	const struct matching_parameters * params,
				unsigned int & id,
				int argc,
				char ** argv )
	{
		if ( obj )
			delete obj;
		obj = new c_recognition;
		matching_process :: params = params;
		if ( obj->setup( 	argc, 
							argv, 
							params->iris_code_fn.c_str(),
							params->fragility_map_fn.c_str(),
							params->nb_fused_images,
							params->nb_thresholds ) )
			return 1;						
		rep = argv[3];
		r_id = &id;
		//Chargement du FBR
		return 0;
	}
	
	int process ( )
	{
		unsigned int i = *r_id;

		if ( i == params->nb_fbr * params->nb_alpha )
			return 1;
		(*r_id) ++;
		double fbr,
			   alpha;
			   
		if ( params->nb_fbr > 1 )
			fbr = ( i % params->nb_fbr ) * ( ( params->fbr_max - params->fbr_min ) / ( params->nb_fbr - 1.0 ) ) + params->fbr_min;
		else
			fbr = params->fbr_min;
		
		if ( params->nb_alpha > 1 )
			alpha = ( i / params->nb_fbr ) * ( ( params->alpha_max - params->alpha_min ) / ( params->nb_alpha - 1.0 ) ) + params->alpha_min;
		else
			alpha = params->alpha_min;	
		
		cout << "Job lauched for " << endl;
		cout << "Distance = " << params->distance << endl;
		cout << "Alpha = " << alpha << endl;
		cout << "FBR = " << fbr << endl;
		cout << "ID = " << i << endl; 
		
		
		if ( obj->match( fbr, alpha ) )
			return 1;
		
		stringstream oss;
		oss << rep << "/" << "Distance=" << params->distance << "&Alpha=" << alpha << "&FBR=" << fbr;
		mkdir( oss.str().c_str(), 014777 );
		obj->save(oss.str().c_str());
		oss.str("");
		return 0;
		
	}
	
	~matching_process()
	{
		if ( obj )
			delete obj;
		obj = NULL;
	}
	
	
	/** @brief
	 * Rep de sauvegarde
	 * 
	 */
	string rep;
	
	/** @brief
	 * Matching
	 * 
	 */
	c_recognition * obj;
	
	/** @brief
	 * Paramètres
	 */
	const struct matching_parameters * params;
	
	/**@brief
	 * Numéro du process.
	 **/
	unsigned int * r_id;
};

void * do_process ( void * data )
{
	matching_process * toto = (matching_process *) data;
	while ( ! toto->process ( ) );
	return NULL;
}


int main( int argc, char ** argv )
{
	matching_parameters m_params;
	api_parameters f_params;
	
	//Chargement des fichiers de paramètres
	for ( int i = 4; i < argc; ++ i )
		f_params.load(argv[i]);
	
	if ( m_params.setup( f_params ) )
		return 1;

	//Construction des différents process
	unsigned int id = 0;
	matching_process * m_process = new matching_process[m_params.nb_process];
	for ( unsigned int i = 0; i < m_params.nb_process; ++ i )
	{
		if ( m_process[i].setup( &m_params, id, argc, argv ) )
		{

			delete[] m_process;
			return 1;
		}

	}

	//Création des threads
	pthread_t * tab_threads = new pthread_t[ m_params.nb_process ];
	for ( unsigned int i = 0; i < m_params.nb_process; ++ i )
	{
		pthread_create( tab_threads + i,
						NULL,
						do_process,
						(void*) ( m_process + i ) ) ;		
	
		usleep( 40000);
	}
	for ( unsigned j = 0; j < m_params.nb_process; ++ j )
		pthread_join( tab_threads[j], NULL);	
	delete[] tab_threads;
	return 0;
		
}
