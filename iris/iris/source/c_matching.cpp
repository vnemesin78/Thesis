#include "c_matching.hpp"

c_matching :: c_matching ( void ) 
: c_database( )
{
	c_matching :: initialize();
}

c_matching :: c_matching ( 	const char * rep_name, 
							const char * iris_code_file,
							const char * fragility_map_file,
							distance::function_prototype_bis dist,
							double d_theta,
							ostream * stream)
: c_database( )
{
	initialize();
	setup ( rep_name, 
			iris_code_file,
			fragility_map_file,
			dist,
			d_theta,
			stream );
}

int c_matching :: setup ( 	const char * rep_name, 
							const char * iris_code_file,
							const char * fragility_map_file,
							distance::function_prototype_bis dist,
							double d_theta,
							ostream * stream)
{
	c_matching :: free();
	c_matching :: initialize();

	c_database :: free();
	c_database :: initialize();
	
	if ( c_database :: setup ( 	rep_name, 
								iris_code_file,
								fragility_map_file,
								stream ) )
		return 1;
	
	//Distance
	c_matching :: dist = dist;
	
	if ( dist == (distance::function_prototype_bis) distance :: Hamming ) 
	{
		dist_params = (void*) new distance::Hamming_parameters;
		dist_bis = distance :: Hamming_opt<unsigned long long>;
	}
	else if ( c_matching :: dist ==  (distance::function_prototype_bis) distance :: fragile_bit_distance )
	{
		dist_params = (void*) new distance::Hamming_parameters;
		dist_bis = distance :: fragile_bit_distance_opt<unsigned long long>;
	}
	else if ( c_matching :: dist == (distance::function_prototype_bis) distance :: Hamming_FBD )
	{
		dist_params = (void*) new distance::Hamming_FBD_parameters;
		dist_bis = distance :: Hamming_FBD_opt<unsigned long long>;
	}
	else if ( c_matching :: dist == (distance::function_prototype_bis) distance :: Hamming_expectation )
	{
		dist_params = NULL;
		dist_bis = NULL;
	}
	else
	{
		if ( err_stream )
			*err_stream << "Error: This distance is not implemented yet!" << endl;
		return 1;
	}
	c_matching :: d_theta = d_theta;
	return 0;
}
	
int c_matching :: setup ( 	const char * rep_name, 
							const char * iris_code_file,
							const char * fragility_map_file,
							api_parameters & params,
							ostream * stream,
							const char * n_space,
							const char * distance_name,
							const char * d_theta_name )	
{
	c_matching :: free();
	c_matching :: initialize();

	c_database :: free();
	c_database :: initialize();
	
	int q = 0;					
	stringstream oss;
	string d;
	distance::function_prototype_bis _dist = 0;
	
	oss << n_space << "::" << distance_name;
	if ( api_get_string( params, 
						 oss.str().c_str(), 
						 &d, 
						 stream ) )
		q = 1;
	oss.str("");
	if ( d == "Hamming" )
		_dist = distance::Hamming;
	else if ( d == "E(Hamming)" )
		_dist = distance::Hamming_expectation;
	else if ( d == "FBD" )
		_dist = distance::fragile_bit_distance;
	else if ( d == "Hamming_FBD" )
		_dist = distance::Hamming_FBD;
	else
	{
		if ( stream )
			*stream << "Error : Unknown distance " << d << endl;
		q = 1;
	}

	double d_theta;
	oss << n_space << "::" << d_theta_name;
	if ( api_get_double( params, 
						 oss.str().c_str(), 
						 &d_theta, 
						 stream ) )
		q = 1;
	oss.str("");
	if ( q )
		return 1;

	return setup( 	rep_name, 
					iris_code_file,
					fragility_map_file,
					_dist,
					d_theta,
					stream  );
}


struct d_matching
{
	unsigned int id;
	double dist;
	bool operator < ( const d_matching & d ) const
	{
		return ( dist < d.dist );
	}	
};

int c_matching :: set_alpha( double alpha )
{
	if ( dist == (distance::function_prototype_bis) distance :: Hamming_FBD )
	{
		( (distance::Hamming_FBD_parameters* ) dist_params)->alpha = alpha;
		return 0;
	}
	return 1;
}


unsigned int c_matching :: matching (	double * distances,
										unsigned int * orders,
										unsigned int & truth,
										const IplImage * iris_code,
										const IplImage * fragility_map,
										const double & fbt,
										const char * class_name,
										const char * video_name )
{
	//Vérité
	truth = _nb_classes;
	for ( unsigned int i = 0; i < _nb_classes; ++ i )
	{
		if ( _class_names[i] == class_name )
		{	
			truth = i;
			break;
		}
	}
	if ( 	dist == (distance::function_prototype_bis) distance :: Hamming 				|| 
			dist == (distance::function_prototype_bis) distance :: fragile_bit_distance 	)
		( (distance::Hamming_parameters* ) dist_params)->fragile_bit_threshold_2 = fbt;
	else if ( dist == (distance::function_prototype_bis) distance :: Hamming_FBD )
		( (distance::Hamming_FBD_parameters* ) dist_params)->p_Hamming.fragile_bit_threshold_2 = fbt;
	 	
	
	//Calcul des distances
	for ( unsigned int i = 0; i < _nb_classes; ++ i)
	{
		if ( 	dist == (distance::function_prototype_bis) distance :: Hamming 				|| 
				dist == (distance::function_prototype_bis) distance :: fragile_bit_distance 	)
			( (distance::Hamming_parameters* ) dist_params)->fragile_bit_threshold_1 = _fragile_bit_thresholds[i];
		else if ( dist == (distance::function_prototype_bis) distance :: Hamming_FBD )
			( (distance::Hamming_FBD_parameters* ) dist_params)->p_Hamming.fragile_bit_threshold_1 = _fragile_bit_thresholds[i];
			
			
		if ( iris_code && _iris_codes[i] )
		{
			//Recadrage
			int theta;
			distance::registering (	distances[i],
									theta,
									_iris_codes[i],
									_fragility_maps[i],
									iris_code,
									fragility_map,
									dist_params,
									dist,
									- (d_theta * _iris_codes[i]->width) / 2,
									+ (d_theta * _iris_codes[i]->width) / 2 );
		}
		else
			distances[i] = 1.0;	
	}
		
	//Tri des distances
	d_matching * tmp = new d_matching[_nb_classes];
	for ( unsigned int i = 0; i < _nb_classes; ++ i )
	{
		tmp[i].id = i;
		tmp[i].dist = distances[i];
	}
	sort ( tmp, tmp + _nb_classes );
	
	for ( unsigned int i = 0; i < _nb_classes; ++ i )
	{
		orders[i] = tmp[i].id;
	}	

	delete[] tmp;
	return orders[0];
}

unsigned int c_matching :: matching_opt (	double * distances,
											unsigned int * orders,
											unsigned int & truth,
											const unsigned long long * iris_code_bis,
											const unsigned long long * mask,
											const unsigned long long * fragile_bits,
											const char * class_name,
											const char * video_name)
{
	

	//Vérité
	truth = _nb_classes;
	for ( unsigned int i = 0; i < _nb_classes; ++ i )
	{
		if ( _class_names[i] == class_name )
		{	
			truth = i;
			break;
		}
	}
	double p = 0;
	if ( 	dist == (distance::function_prototype_bis) distance :: Hamming_FBD 	)
	{
		p = ( (distance::Hamming_FBD_parameters* ) dist_params)->alpha;
	}
	
	//Calcul des distances
	
	unsigned long long * _code_2_data = new unsigned long long[width_step * 40];
	unsigned long long * _mask_2_data = new unsigned long long[width_step * 40];
	unsigned long long * _fb_2_data = new unsigned long long[width_step * 40];
			
	
	for ( unsigned int i = 0; i < _nb_classes; ++ i)
	{
		if (  iris_code_bis &&  _iris_codes_bis[i] )
		{

			
			//Recadrage
			int theta;
			distance::registering_64b(	distances[i],
										theta,
										_iris_codes_bis[i],
										_masks[i],
										_fragile_bits[i],
										iris_code_bis,
										mask,
										fragile_bits,
										_iris_codes[i]->width,
										_iris_codes[i]->height,
										width_step,
										- (d_theta * _iris_codes[i]->width) / 2,
										+ (d_theta * _iris_codes[i]->width) / 2,
										dist_bis,
										(const void*) &p,
										_code_2_data,
										_mask_2_data,
										_fb_2_data
								 );	
		}
		else
			distances[i] = 1.0;	
	}
	delete[] _code_2_data;
	delete[] _mask_2_data;
	delete[] _fb_2_data;
		
	//Tri des distances
	d_matching * tmp = new d_matching[_nb_classes];
	for ( unsigned int i = 0; i < _nb_classes; ++ i )
	{
		tmp[i].id = i;
		tmp[i].dist = distances[i];
	}
	sort ( tmp, tmp + _nb_classes );
	
	for ( unsigned int i = 0; i < _nb_classes; ++ i )
	{
		orders[i] = tmp[i].id;
	}	

	delete[] tmp;
	return orders[0];
}






c_matching :: ~c_matching()
{
	free();
	initialize();
}

void c_matching :: initialize()
{
	dist = 0;
	dist_params = 0;	
	d_theta = 0;	
	//Distances

}

void c_matching :: free()
{
	if ( dist_params )
	{
		if ( dist == (distance::function_prototype_bis) distance :: Hamming ) 
			delete (distance::Hamming_parameters*) dist_params;
		else if ( c_matching :: dist ==  (distance::function_prototype_bis) distance :: fragile_bit_distance )
			delete (distance::Hamming_parameters*) dist_params;
		else if ( c_matching :: dist == (distance::function_prototype_bis) distance :: Hamming_FBD )
			 delete (distance::Hamming_FBD_parameters*) dist_params;

	}
}

