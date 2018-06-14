#include "c_fusion_iris_template.hpp"
#include <sys/stat.h>
c_fusion_iris_template :: c_fusion_iris_template(  void	)
{
	initialize();
}

c_fusion_iris_template :: c_fusion_iris_template( 	unsigned int nb_images,
													int type,
													unsigned int nb_directions,
													unsigned int nb_samples,
													double d_theta,
													ostream * stream )
{
	initialize();
	setup ( nb_images, 
			type,
			nb_directions, 
			nb_samples, 
			d_theta,
			stream );
}

int c_fusion_iris_template :: setup ( 	unsigned int nb_images,
										int type,
										unsigned int nb_directions,
										unsigned int nb_samples,
										double d_theta,
										ostream * stream )
{
	free();
	initialize();
	err_stream = stream;
	
	//Check des arguments
	if (	! nb_images 	|| 
			!nb_directions 	|| 
			!nb_samples 	)
	{
		if ( stream )
			*stream << "Error in int c_fusion_iris_template :: setup : ( ! nb_images || !nb_directions || !nb_samples ) " << endl;
		return 1;
	}
	
	_nb_images_max = nb_images;
	_type = type;
	_nb_directions = nb_directions;
	_nb_samples = nb_samples;
	c_fusion_iris_template :: d_theta = d_theta;
	alloc();
	memset( _scores, 
			0, 
			sizeof( double ) * _nb_images_max );
	_fusion_obj = new c_fusion;
	if ( _fusion_obj->setup ( 	nb_directions, 
								2 * nb_samples, 
								stream ) )
		return 1;
	return 0;
}

void c_fusion_iris_template :: alloc()
{
	_iris_codes = new IplImage*[_nb_images_max];
	_iris_code_masks = new IplImage*[_nb_images_max];
	_scores = new double[_nb_images_max];
	for ( unsigned int i = 0; i < _nb_images_max; ++ i )
	{
		_iris_codes[i] = 
			cvCreateImage( 	cvSize( _nb_directions, 
									_nb_samples * 2 ),
							IPL_DEPTH_8U,
							1 );
		_iris_code_masks[i] = 
			cvCreateImage( 	cvSize( _nb_directions, 
									_nb_samples * 2 ),
							IPL_DEPTH_8U,
							1 );				
	}
	
	_saved_fragility_maps = new IplImage*[_nb_images_max];
	_saved_codes = new IplImage*[_nb_images_max];
	for ( unsigned int i = 0; i < _nb_images_max; ++ i )
	{
		_saved_fragility_maps[i] = 
			cvCreateImage( 	cvSize( _nb_directions, 
									_nb_samples * 2 ),
							IPL_DEPTH_8U,
							1 );
		_saved_codes[i] = 
			cvCreateImage( 	cvSize( _nb_directions, 
									_nb_samples * 2 ),
							IPL_DEPTH_8U,
							1 );	
		
		
		
	}
	
}

int c_fusion_iris_template :: setup ( 	api_parameters & params,
										const char * n_space,
										const char * nb_images_n,
										const char * type_n,
										const char * n_space_iris_code,
										const char * nb_dir_name,
										const char * nb_samples_name,	
										const char * d_theta_name,				
										ostream * stream )
{
	free();
	initialize();
	
	err_stream = stream;	
	stringstream oss;
	int q = 0;
	
	
	oss << n_space << "::" << nb_images_n;
	if ( api_get_positive_integer( 	params, 
									oss.str().c_str(), 
									&_nb_images_max, 
									stream ) )
		q = 1;
	oss.str("");
	
	oss << n_space << "::" << type_n;
	if ( api_get_integer( 	params, 
							oss.str().c_str(), 
							&_type, 
							stream ) )
		q = 1;
	oss.str("");	
	
	
	if ( q )
		return 1;
	

		
	oss << n_space_iris_code << "::" << nb_dir_name;
	if ( api_get_positive_integer( 	params, 
									oss.str().c_str(), 
									&_nb_directions, 
									stream ) )
		q = 1;
	oss.str("");	
	
	oss << n_space_iris_code << "::" << nb_samples_name;
	if ( api_get_positive_integer( 	params, 
									oss.str().c_str(), 
									&_nb_samples, 
									stream ) )
		q = 1;
	oss.str("");	
		
	double d_theta;
	oss << n_space_iris_code << "::" << d_theta_name;
	if ( api_get_double( 	params, 
							oss.str().c_str(), 
							&d_theta, 
							stream ) )
		q = 1;
	oss.str("");	
		
		
	if ( q )
		return 1;
	return  setup( 	_nb_images_max, 
					_type,
					_nb_directions, 
					_nb_samples, 
					d_theta,
					stream );	
}


c_fusion_iris_template :: ~c_fusion_iris_template()
{
	free();
	initialize();
}





void c_fusion_iris_template :: initialize()
{
	err_stream = 0;
	_type = 0;
	_nb_directions = 0;
	_nb_samples = 0;
	_nb_images = 0;
	_nb_images_max = 0;
	_iris_codes = 0;
	_iris_code_masks = 0;
	_scores = 0;
	_fusion_obj = 0;
	d_theta = 0;
	_saved_codes = 0;
	_saved_fragility_maps = 0;
}

void c_fusion_iris_template :: free()
{
	
	if ( _iris_codes )
	{
		for ( unsigned int i = 0; i < _nb_images_max; ++ i )
		{
			cvReleaseImage( _iris_codes + i );
		}
		delete[] _iris_codes;
	}
	
	if ( _iris_code_masks )
	{
		for ( unsigned int i = 0; i < _nb_images_max; ++ i )
		{
			cvReleaseImage( _iris_code_masks + i );
		}
		delete[] _iris_code_masks;
	}
		
	if ( _saved_codes )
	{
		for ( unsigned int i = 0; i < _nb_images_max; ++ i )
		{
			cvReleaseImage( _saved_codes + i );
		}
		delete[] _saved_codes;
	}
	
	if ( _saved_fragility_maps )
	{
		for ( unsigned int i = 0; i < _nb_images_max; ++ i )
		{
			cvReleaseImage( _saved_fragility_maps + i );
		}
		delete[] _saved_fragility_maps;
	}
	
	if ( _scores )
		delete[] _scores;
	
	if ( _fusion_obj )
		delete _fusion_obj;
	
}

void c_fusion_iris_template :: add_iris( 	const IplImage * code,
											const IplImage * mask,
											double score )
{
	//Données tampons
	IplImage * p_code, 
			  *p_mask;
	for ( unsigned int i = 0; i < _nb_images_max; ++ i )
	{
		if ( _nb_images < _nb_images_max )
			_nb_images ++;
		
		
		if ( score > _scores[i] )
		{
			p_code = _iris_codes[_nb_images_max - 1];
			p_mask = _iris_code_masks[_nb_images_max - 1];
			
			for ( unsigned int j = _nb_images_max - 1; j > i; --j )
			{
				_iris_codes[j] = _iris_codes[ j - 1 ];
				_iris_code_masks[j] = _iris_code_masks[ j - 1 ];
				_scores[j] = _scores[ j -  1 ];
			}  
			_iris_codes[i] = p_code;
			cvCopyImage ( code, _iris_codes[i] );
			
			_iris_code_masks[i] = p_mask;
			cvSetImageROI( _iris_code_masks[i],
						   cvRect( 0,0, mask->width, mask->height ) );
			cvCopyImage ( mask, _iris_code_masks[i] );		
			cvSetImageROI( _iris_code_masks[i],
						   cvRect( 0,mask->height, mask->width, mask->height ) );
			cvCopyImage ( mask, _iris_code_masks[i] );	
			cvResetImageROI(_iris_code_masks[i]);
			
			_scores[i] = score;
			break;
		}
	}
	
}

void c_fusion_iris_template :: get_iris_code_and_fragility_map( IplImage * code, 
																IplImage * map )
{
	for (unsigned int row = 0; row < _nb_samples; ++ row )
	{
		for (unsigned int col = 0; col < _nb_directions; ++ col )
		{
			if ( GET_IMAGE_PIXEL(	map, 
									row, 
									col, 
									unsigned char ) )
			{
				GET_IMAGE_PIXEL( map, row, col, unsigned char ) 
					= abs( 255 - 2 * ( (int) GET_IMAGE_PIXEL(code, row, col, unsigned char ) ) );
				GET_IMAGE_PIXEL( map, row + _nb_samples, col, unsigned char ) 
					= abs( 255 - 2 * ( (int) GET_IMAGE_PIXEL(code, row + _nb_samples, col, unsigned char ) ) );
				if ( GET_IMAGE_PIXEL(	code, row, col, unsigned char ) > 127 )
					GET_IMAGE_PIXEL( code, row, col, unsigned char ) = 255;
				else
					GET_IMAGE_PIXEL( code, row, col, unsigned char ) = 0;
					
				if ( GET_IMAGE_PIXEL(	code, row + _nb_samples, col, unsigned char ) > 127 )
					GET_IMAGE_PIXEL( code, row + _nb_samples, col, unsigned char ) = 255;
				else
					GET_IMAGE_PIXEL( code, row + _nb_samples, col, unsigned char ) = 0;
					
			}
			else
			{
				GET_IMAGE_PIXEL(map, row, col, unsigned char ) = 0;
				GET_IMAGE_PIXEL(map, row + _nb_samples, col, unsigned char ) = 0;
				GET_IMAGE_PIXEL( code, row + _nb_samples, col, unsigned char ) = 128;
				GET_IMAGE_PIXEL( code, row, col, unsigned char ) = 128;
			}
		}	
	}
	
	
	
}



void c_fusion_iris_template :: fusion_data ( void )
{
	distance::Hamming_parameters dist_params;
	dist_params.fragile_bit_threshold_1 = 1;
	dist_params.fragile_bit_threshold_2 = 1;
	
	double score;
	IplImage * _fuzzy_code = _saved_codes[0];
	IplImage * _fragility_map = _saved_fragility_maps[0];
	
	//0
	_fusion_obj->reset();
	score = _scores[0];
	switch ( _type )
	{
		case(1): //exp mean
			score = exp( 1 );
		break;
		case(2): //mean
			score = 1;
		break;
		case(3):
			//Mauvaise images
			if ( score < 0.2 ) 
				score = 0.5 * score;
			//Images moyennes
			else if ( score >= 0.2 && score < 0.4 )
				score = 0.1 + 0.8 * (score - 0.2) / 0.2;
			else
				score = 0.9 + 0.1 * (score - 0.4) / 0.6;
		break;
		default://weighted mean
			score = 1;
		break;
		
	}		

	_fusion_obj->add_image( _iris_codes[0], 
							_iris_code_masks[0],
							score,
							0 );	
	_fusion_obj->compute (	_fuzzy_code,
							_fragility_map,
							true,
							0.5 );
							
	//Calcul de la carte de fragilité
	get_iris_code_and_fragility_map( 	_fuzzy_code, 
										_fragility_map );

	
	
	for ( unsigned int i = 1; i < _nb_images; ++ i )
	{
		double d = 0;
		int theta = 0;

		//Registering
		distance:: registering (	d,
									theta,
									_fuzzy_code,
									_fragility_map,
									_iris_codes[i],
									_iris_code_masks[i],
									(void*) &dist_params,
									distance::Hamming,
									- (d_theta * _iris_codes[i]->width) / 2,
									(d_theta * _iris_codes[i]->width) / 2 );	
									
		//~ distance:: registering (	d,
									//~ theta,
									//~ _saved_codes[0],
									//~ _saved_fragility_maps[0],
									//~ _iris_codes[i],
									//~ _iris_code_masks[i],
									//~ (void*) &dist_params,
									//~ distance::Hamming,
									//~ - (d_theta * _iris_codes[i]->width) / 2,
									//~ (d_theta * _iris_codes[i]->width) / 2 );	
												
		//Sauvegarde				
		_fuzzy_code = _saved_codes[i];
		_fragility_map = _saved_fragility_maps[i];


		score = _scores[i];
		switch ( _type )
		{
			case(1): //exp mean
				score = _scores[i] / _scores[0];
				score = exp( score );
			break;
			case(2): //mean
				score = 1;
			break;
			case(3):
				//Mauvaise images
				if ( score < 0.2 ) 
					score = 0.5 * score;
				//Images moyennes
				else if ( score >= 0.2 && score < 0.4 )
					score = 0.1 + 0.8 * (score - 0.2) / 0.2;
				else
					score = 0.9 + 0.1 * (score - 0.4) / 0.6;
			break;
			
			
			default://weighted mean
			
			
			
			break;
			
		}		
		
		//Ajout de l'image
		_fusion_obj->add_image( _iris_codes[i],
								_iris_code_masks[i],
								score,
								theta );
								
	
		//Calcul de la fusion à la fin
		_fusion_obj->compute (	_fuzzy_code,
								_fragility_map,
								true,
								0.5 );
		get_iris_code_and_fragility_map( 	_fuzzy_code, 
											_fragility_map );
	}
}

int c_fusion_iris_template :: fusion ( 	const iris_data * i_data,
											unsigned int nb_iris )
{
	memset ( _scores, 0, sizeof(double) * _nb_images_max );
	_nb_images = 0;
	for ( unsigned int i = 0; i < nb_iris; ++ i )
	{
		add_iris ( 	i_data[i].code,
					i_data[i].code_mask,
					i_data[i].nrj_ratio );
		
	}
	
	if ( _nb_images == 0 )
	{
		if ( err_stream )
			*err_stream << "No iris code for the class!" << endl;
		return 1;
	}
	fusion_data();
	return 0;
}

int c_fusion_iris_template :: fusion ( const char * rep )
{
	memset ( _scores, 0, sizeof(double) * _nb_images_max );
	_nb_images = 0;
	FILE * file = NULL;
	stringstream oss;
	unsigned int nb = 0;	
	do
	{
		oss.str("");
		iris_data i_data;

		if ( file )
			fclose(file);
		oss << rep << "/iris_code_" << nb <<".png";
		if  ( nb >= 1 )
		{
			if ( ! i_data.load( rep, nb - 1 ) )
			{
				add_iris ( 	i_data.code,
							i_data.code_mask,
							i_data.nrj_ratio );
				
			}
		}
		nb ++;
	//~ } while ( ( file = fopen(oss.str().c_str(), "r" ) ) != NULL || );
		file = fopen(oss.str().c_str(), "r" );
	} while ( nb < 100 );
	if ( _nb_images == 0 )
	{
		if ( err_stream )
			*err_stream << "No iris code for the class!" << endl;
		return 1;
	}
	fusion_data();	
	return 0;
}
						
int c_fusion_iris_template :: save ( 	const char * rep,
										const char * code_n,
										const char * fragil_n  )
{
#ifdef __linux
	mkdir ( rep, 014777 );
#elif __APPLE__ 
	mkdir ( rep, 014777 );
#elif _WIN32
	_mkdir( rep );
#elif _WIN64
	_mkdir( rep );
#else
	#error
#endif
	if ( _nb_images == 0 )
		return 1;
	
	stringstream oss;
	
	
	for ( unsigned int i = 0; i < _nb_images; ++ i )
	{
		oss << rep << "/" << code_n << "_" << i + 1 << ".png";
		cvSaveImage( oss.str().c_str(), _saved_codes[i] );
		oss.str("");
		
		oss << rep << "/" << fragil_n << "_" << i + 1 << ".png";
		cvSaveImage( oss.str().c_str(), _saved_fragility_maps[i] );
		oss.str("");
		
	}
	for ( unsigned int i =  _nb_images; i < _nb_images_max; ++ i )
	{
		oss << rep << "/" << code_n << "_" << i + 1 << ".png";
		cvSaveImage( oss.str().c_str(), _saved_codes[ _nb_images - 1] );
		oss.str("");
		
		oss << rep << "/" << fragil_n << "_" << i + 1 << ".png";
		cvSaveImage( oss.str().c_str(), _saved_fragility_maps[ _nb_images - 1] );
		oss.str("");
	}
	
	
	return 0;
}
