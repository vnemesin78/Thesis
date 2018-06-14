#include "mouloud.hpp"
#include "image_utility.hpp"	
#include <cstdlib>
#include <opencv/highgui.h>
Mouloud :: Mouloud ( void )
{
	initialize();
}

Mouloud :: Mouloud (	unsigned int nb_iter_gem,
						unsigned int nb_iter_icm,
						pdf_type bg_function,
						estimation_type est_bg_function,		
						const void * est_bg_params,
						unsigned int size_est_bg_params,
						pdf_type obj_function,
						estimation_type est_obj_function,
						const void * est_obj_params,
						unsigned int size_est_obj_params,			
						pdf_field_type edge_function,
						const void * edge_params,
						unsigned int size_edge_params,
						unsigned int width,
						unsigned int height)
{
	initialize();
	setup (	nb_iter_gem,
			nb_iter_icm,
			bg_function,
			est_bg_function,
			est_bg_params,
			size_est_bg_params,
			obj_function,
			est_obj_function,
			est_obj_params,
			size_est_obj_params,
			edge_function,
			edge_params,
			size_edge_params,
			width,
			height );
	
	
	
	
	
}

int Mouloud :: setup ( )
{
	_free();
	initialize();
	return 0;
}

int Mouloud :: setup (	unsigned int nb_iter_gem,
						unsigned int nb_iter_icm,
						pdf_type bg_function,
						estimation_type est_bg_function,		
						const void * est_bg_params,
						unsigned int size_est_bg_params,
						pdf_type obj_function,
						estimation_type est_obj_function,
						const void * est_obj_params,
						unsigned int size_est_obj_params,			
						pdf_field_type edge_function,
						const void * edge_params,
						unsigned int size_edge_params,
						unsigned int width,
						unsigned int height )
{
	_free();
	initialize();
	
	_nb_iter_gem = nb_iter_gem;
	_nb_iter_icm = nb_iter_icm;
	
	_bg_function = bg_function;
	_est_bg_function = est_bg_function;
	
	_size_est_bg_params = size_est_bg_params;
	
	_obj_function = obj_function;
	_est_obj_function = est_obj_function;
	
	_size_est_obj_params = size_est_obj_params;	
	
	_edge_function = edge_function;
	_size_edge_params = size_edge_params;
	
	_width = width;
	_height = height;
	
	alloc();
	memcpy (	_edge_params, 
				edge_params, 
				_size_edge_params );
	memcpy ( 	_est_obj_params,
				est_obj_params,
				_size_est_obj_params );
	memcpy (	_est_bg_params,
				est_bg_params,
				_size_est_bg_params );
	return 0;
}

Mouloud :: ~Mouloud()
{
	_free();
	initialize();
}

void Mouloud :: _free()
{
	if ( _log_bg_pdf )
		cvReleaseImage( &_log_bg_pdf );
	if ( _log_obj_pdf )
		cvReleaseImage( &_log_obj_pdf );		
	if ( _log_r_pdf )
		cvReleaseImage( &_log_r_pdf );
	if ( _mask )
		cvReleaseImage( &_mask );	
	if ( _tmp_image )
		cvReleaseImage ( &_tmp_image );
	free( _est_bg_params );
	free( _est_obj_params );
	free( _edge_params );
}

void Mouloud :: initialize()
{
	err_stream = &cout;

	_width = 0;
	_height = 0;

	_log_bg_pdf = 0;
	_log_obj_pdf = 0;
	_log_r_pdf = 0;
	_tmp_image = 0;
	_mask = 0;

	//Params
	_nb_iter_icm = 0;
	_nb_iter_gem = 0;

	_bg_function = 0;
	_est_bg_function = 0;		
	_est_bg_params = 0;
	_size_est_bg_params = 0;
	
	_obj_function = 0;
	_est_obj_function = 0;
	_est_obj_params = 0;
	_size_est_obj_params = 0;
		
	//Edge
	_edge_function = 0;
	_edge_params = 0;
	_size_edge_params = 0;
}

void Mouloud :: alloc()
{
	if ( _width && _height )
	{
		_log_bg_pdf = cvCreateImage ( cvSize( _width, _height + 1 ),
									  IPL_DEPTH_64F,
									  1 );
		_log_obj_pdf = cvCreateImage (	cvSize( _width, _height + 1 ),
										IPL_DEPTH_64F,
										1 );
		_log_r_pdf = cvCreateImage (	cvSize( _width, _height + 1 ),
										IPL_DEPTH_64F,
										1 );
		_mask = cvCreateImage (	cvSize( _width, _height ),
								IPL_DEPTH_8U,
								1 );
		_tmp_image = cvCreateImage (	cvSize( _width, _height ),
										IPL_DEPTH_64F,
										1);
	}
	
	if ( _size_est_bg_params )
		_est_bg_params = malloc( _size_est_bg_params );
	if ( _size_est_obj_params )
		_est_obj_params = malloc ( _size_est_obj_params );
	if ( _size_edge_params )
		_edge_params = malloc ( _size_edge_params );
	
	
	
	
	
	
	
	
	
	
	
}

int Mouloud :: compute_edges( 	double * radii,
								void * bg_params,
								void * obj_params,
								const IplImage * image,
								const IplImage * mask )
{
	unsigned int 	x_offset,
					y_offset,
					width,
					height,
					width_step;
	
	unsigned int 	m_x_offset,
					m_y_offset,
					m_width,
					m_height,
					m_width_step;
	
	
	
	if ( ! radii || ! bg_params || !obj_params || !image || !mask )
	{
		*err_stream << "Error : Invalid argument(s) in int Mouloud :: compute_edges" << endl;
		return -1;
	}

	GET_IMAGE_DIM( 	image, 
					width, 
					height, 
					width_step,
					x_offset, 
					y_offset );
	
	GET_IMAGE_DIM( 	mask, 
					m_width, 
					m_height, 
					m_width_step,
					m_x_offset, 
					m_y_offset );
	
	
	if ( 	mask->depth != IPL_DEPTH_8U 		||
			m_width != width					||
			m_height != height 					)
	{
		*err_stream << "Error : Invalid mask in int Mouloud :: compute_edges" << endl;
		return -1;
	}	
	switch ( image->depth )	
	{
		case( IPL_DEPTH_8S):
			return ( compute_edges(	radii,
										bg_params,
										obj_params,
										( (char *) ( image->imageData + y_offset * width_step ) ) + x_offset,
										( (unsigned char *) ( mask->imageData + m_y_offset * m_width_step ) ) + m_x_offset,
										width,
										height,
										width_step / sizeof( char ),
										m_width_step / sizeof( unsigned char ) )
						);
		break;
		case( IPL_DEPTH_8U):
			return ( compute_edges(	radii,
										bg_params,
										obj_params,
										( (unsigned char *) ( image->imageData + y_offset * width_step ) ) + x_offset,
										( (unsigned char *) ( mask->imageData + m_y_offset * m_width_step ) ) + m_x_offset,
										width,
										height,
										width_step / sizeof( unsigned char ),
										m_width_step / sizeof( unsigned char ) )
						);
		
		
		
		break;
		case( IPL_DEPTH_16S):
			return ( compute_edges(	radii,
										bg_params,
										obj_params,
										( (short int *) ( image->imageData + y_offset * width_step ) ) + x_offset,
										( (unsigned char *) ( mask->imageData + m_y_offset * m_width_step ) ) + m_x_offset,
										width,
										height,
										width_step / sizeof( short int ),
										m_width_step / sizeof( unsigned char ) )
						);
		
		
		
		break;		
		case( IPL_DEPTH_16U):
			return ( compute_edges(	radii,
										bg_params,
										obj_params,
										( (unsigned short int *) ( image->imageData + y_offset * width_step ) ) + x_offset,
										( (unsigned char *) ( mask->imageData + m_y_offset * m_width_step ) ) + m_x_offset,
										width,
										height,
										width_step / sizeof( unsigned short int ),
										m_width_step / sizeof( unsigned char ) )
						);		
		break;		
		case( IPL_DEPTH_32S):
			return ( compute_edges(	radii,
										bg_params,
										obj_params,
										( (long int *) ( image->imageData + y_offset * width_step ) ) + x_offset,
										( (unsigned char *) ( mask->imageData + m_y_offset * m_width_step ) ) + m_x_offset,
										width,
										height,
										width_step / sizeof( long int ),
										m_width_step / sizeof( unsigned char ) )
						);		
		break;		
		case( IPL_DEPTH_32F):

			return ( compute_edges(	radii,
										bg_params,
										obj_params,
										( (float *) ( image->imageData + y_offset * width_step ) ) + x_offset,
										( (unsigned char *) ( mask->imageData + m_y_offset * m_width_step ) ) + m_x_offset,
										width,
										height,
										width_step / sizeof( float ),
										m_width_step / sizeof( unsigned char ) )
						);		
		break;		
		case( IPL_DEPTH_64F):

			return ( compute_edges(	radii,
										bg_params,
										obj_params,
										( (double *) ( image->imageData + y_offset * width_step ) ) + x_offset,
										( (unsigned char *) ( mask->imageData + m_y_offset * m_width_step ) ) + m_x_offset,
										width,
										height,
										width_step / sizeof( double ),
										m_width_step / sizeof( unsigned char ) )
						);		
		break;		
		default:
			*err_stream << "Error : Unsupported image format in int Mouloud :: compute_edges" << endl;
			return -1;
		break;
	}
}

template<class type> 
int Mouloud :: compute_edges( 	double * radii,
								void * bg_params,
								void * obj_params,
								const type * data,
								const unsigned char * data_mask,
								unsigned int width,
								unsigned int height,
								unsigned int width_step,
								unsigned int mask_width_step )
{
	if ( ! radii || ! bg_params || !obj_params || !data || !data_mask || !width || !height )
	{
		*err_stream << "Error : Invalid argument(s) in int Mouloud :: compute_edges" << endl;
		return -1;
	}

	
	if ( width >_width or height > _height )
	{
		if ( _log_bg_pdf )
			cvReleaseImage( &_log_bg_pdf );
		if ( _log_obj_pdf )
			cvReleaseImage( &_log_obj_pdf );		
		if ( _log_r_pdf )
			cvReleaseImage( &_log_r_pdf );
		if ( _mask )
			cvReleaseImage( &_mask );	
		if ( _tmp_image )
			cvReleaseImage (&_tmp_image );
			
		_width = width;
		_height = height;
		
		if ( _width && _height )
		{
			_log_bg_pdf = cvCreateImage ( cvSize( _width, _height + 1 ),
										  IPL_DEPTH_64F,
										  1 );
			_log_obj_pdf = cvCreateImage (	cvSize( _width, _height + 1 ),
											IPL_DEPTH_64F,
											1 );
			_log_r_pdf = cvCreateImage (	cvSize( _width, _height + 1 ),
											IPL_DEPTH_64F,
											1 );
			_mask = cvCreateImage (	cvSize( _width, _height ),
									IPL_DEPTH_8U,
									1 );
			_tmp_image = cvCreateImage (	cvSize( _width, _height ),
											IPL_DEPTH_64F,
											1 );
		}
			
	} 
		
	for ( unsigned int i = 0; i < height; ++ i )
	{
		for ( unsigned int j = 0; j < width; ++ j )
		{
			( (double*) ( _tmp_image->imageData + i * _tmp_image->widthStep ) )[j] 
				= (double) data[ i * width_step + j ]; 
		}
	}
	return compute_edges( 	radii, 
							bg_params, 
							obj_params,  
							(double*) ( _tmp_image->imageData ), 
							data_mask,
							width,
							height,
							_tmp_image->widthStep / sizeof(double),
							mask_width_step );
}

int Mouloud :: compute_edges( 	double * radii,
								void * bg_params,
								void * obj_params,
								const double * data,
								const unsigned char * data_mask,
								unsigned int width,
								unsigned int height,
								unsigned int width_step,
								unsigned int mask_width_step )
{
	if ( ! radii || ! bg_params || !obj_params || !data || !data_mask || !width || !height )
	{
		*err_stream << "Error : Invalid argument(s) in int Mouloud :: compute_edges" << endl;
		return -1;
	}

	if ( width >_width or height > _height )
	{
		if ( _log_bg_pdf )
			cvReleaseImage( &_log_bg_pdf );
		if ( _log_obj_pdf )
			cvReleaseImage( &_log_obj_pdf );		
		if ( _log_r_pdf )
			cvReleaseImage( &_log_r_pdf );
		if ( _mask )
			cvReleaseImage( &_mask );	
		if ( _tmp_image )
			cvReleaseImage (&_tmp_image );
			
		_width = width;
		_height = height;
		
		if ( _width && _height )
		{
			_log_bg_pdf = cvCreateImage ( cvSize( _width, _height + 1 ),
										  IPL_DEPTH_64F,
										  1 );
			_log_obj_pdf = cvCreateImage (	cvSize( _width, _height + 1 ),
											IPL_DEPTH_64F,
											1 );
			_log_r_pdf = cvCreateImage (	cvSize( _width, _height + 1 ),
											IPL_DEPTH_64F,
											1 );
			_mask = cvCreateImage (	cvSize( _width, _height ),
									IPL_DEPTH_8U,
									1 );
			_tmp_image = cvCreateImage (	cvSize( _width, _height ),
											IPL_DEPTH_64F,
											1 );			
		
		}
			
	} 
		
		
		
		
		
		
		
		
		
		
	for ( unsigned int i = 0; i < _nb_iter_gem; ++ i )
	{
		//Calcul des nouveaux paramètres
		if ( estimate_pdf (	radii,
							bg_params,
							obj_params,
							data,
							data_mask,
							width,
							height,
							width_step,
							mask_width_step ) )
		{
			//~ *err_stream << "Warning : Edges computation failed!" << endl;
			return 1;
		}
		//ICM
		
		if ( compute_bg_pdf(	radii,
								bg_params,
								obj_params,
								data,
								data_mask,
								width,
								height,
								width_step,
								mask_width_step ) )
		{
			//~ *err_stream << "Warning : Edges computation failed!" << endl;
			return 1;
		}
		if ( compute_obj_pdf(	radii,
								bg_params,
								obj_params,
								data,
								data_mask,
								width,
								height,
								width_step,
								mask_width_step ) )
		{
			//~ *err_stream << "Warning : Edges computation failed!" << endl;
			return 1;
		}
	
		if ( compute_radii_pdf(	radii,
								bg_params,
								obj_params,
								data,
								data_mask,
								width,
								height,
								width_step,
								mask_width_step ) )
		{
			//~ *err_stream << "Warning : Edges computation failed!" << endl;
			return 1;
		}	
	
		for ( unsigned int j = 0; j < _nb_iter_icm; ++ j )
		{

								
								
								
										
			if ( estimate_radii(	radii,
									bg_params,
									obj_params,
									data,
									data_mask,
									width,
									height,
									width_step,
									mask_width_step ) )
			{
				//~ *err_stream << "Warning : Edges computation failed!" << endl;
				return 1;
			}
		}
	}
	
	return 0;
}

int Mouloud :: compute_radii_pdf( 	double * radii,
									void * bg_params,
									void * obj_params,
									const double * data,
									const unsigned char * data_mask,
									unsigned int width,
									unsigned int height,
									unsigned int width_step,
									unsigned int mask_width_step )
{
	for ( unsigned int i = 0; i <= height; ++ i )
	{
		for ( unsigned int j = 0; j < width; ++ j )
		{
			double _r = radii[ (j - 1) % width ],
					r = i,
					r_ = radii[ (j + 1) % width ];
			( (double*) (_log_r_pdf->imageData + _log_r_pdf->widthStep * i ) )[j] 
				= _edge_function( _r, r, r_, _edge_params ) 
				+ ( (double*) (_log_bg_pdf->imageData + _log_bg_pdf->widthStep * i ) )[j]
				+ ( (double*) (_log_obj_pdf->imageData + _log_obj_pdf->widthStep * i ) )[j];
		}
	}
	return 0;
}
	
int Mouloud :: compute_bg_pdf(	double * radii,
								void * bg_params,
								void * obj_params,
								const double * data,
								const unsigned char * data_mask,
								unsigned int width,
								unsigned int height,
								unsigned int width_step,
								unsigned int mask_width_step )
{
	//Mise à zéro de la première ligne
	memset(	_log_bg_pdf->imageData + _log_bg_pdf->widthStep * height,
			0,
			_log_bg_pdf->widthStep );
	for ( unsigned int i = 1; i <= height; ++ i )
	{
		unsigned int r = height - i;
		for ( unsigned int j = 0; j < width; ++ j )
		{
			( ( double *) ( _log_bg_pdf->imageData + _log_bg_pdf->widthStep * r ) )[j]
				= ( ( double *) ( _log_bg_pdf->imageData + _log_bg_pdf->widthStep * ( r + 1 ) ) )[j];
			
			if (  data_mask [ r * mask_width_step + j ] )
			{
				double v = data[ r * width_step + j ];
				( ( double *) ( _log_bg_pdf->imageData + _log_bg_pdf->widthStep * r ) )[j] 
					+= _bg_function (	v, bg_params );
			}
			else
				( ( double *) ( _log_bg_pdf->imageData + _log_bg_pdf->widthStep * r ) )[j]  += log(2);
		}
	}
	
	return 0;
}
	
int Mouloud :: compute_obj_pdf(	double * radii,
								void * bg_params,
								void * obj_params,
								const double * data,
								const unsigned char * data_mask,
								unsigned int width,
								unsigned int height,
								unsigned int width_step,
								unsigned int mask_width_step )
{
	//Mise à zéro de la première ligne
	memset(	_log_obj_pdf->imageData,
			0,
			_log_obj_pdf->widthStep );
	for ( unsigned int r = 1; r <= height; ++ r )
	{
		for ( unsigned int j = 0; j < width; ++ j )
		{
			( ( double *) ( _log_obj_pdf->imageData + _log_obj_pdf->widthStep * r ) )[j]
				= ( ( double *) ( _log_obj_pdf->imageData + _log_obj_pdf->widthStep * ( r - 1 ) ) )[j];
			
			if (  data_mask [ ( r - 1 ) * mask_width_step + j ] )
			{
				double v = data[ ( r - 1 ) * width_step + j ];
				( ( double *) ( _log_obj_pdf->imageData + _log_obj_pdf->widthStep * r ) )[j] 
					+= _obj_function (	v, obj_params );
			}
			else
				( ( double *) ( _log_bg_pdf->imageData + _log_bg_pdf->widthStep * r ) )[j]  += log(2);
		}
	}
	
	return 0;
}	
	
	
int Mouloud :: estimate_pdf( 	double * radii,
								void * bg_params,
								void * obj_params,
								const double * data,
								const unsigned char * data_mask,
								unsigned int width,
								unsigned int height,
								unsigned int width_step,
								unsigned int mask_width_step )
{
	//Construction des masques
	//1 pour l'obj.
	//2 pour le BG
	//0 valeur incorrecte
	
	memset( _mask->imageData,
			0,
			_mask->widthStep * height );
			
	for ( unsigned int i = 0; i < height; ++ i )
	{
		for ( unsigned int j = 0; j < width; ++ j )
		{
			if ( data_mask[i * mask_width_step + j] == 255 )
			{
				if ( i < radii[j] )
				( (unsigned char*) (_mask->imageData + i * _mask->widthStep) )[j]
					= 128;
				else
				( (unsigned char*) (_mask->imageData + i * _mask->widthStep) )[j]
					= 254;
			}
		}
	}

	//Estimation
	if ( _est_obj_function(	obj_params,
							data,
							(unsigned char *) _mask->imageData,
							width,
							height,
							width_step,
							_mask->widthStep,
							128,
							_est_obj_params ) )
	{
		//~ *err_stream << "Warning : Obj. pdf estimation!" << endl;
		return 1;
	}

	
	if ( _est_bg_function(	bg_params,
							data,
							(unsigned char *) _mask->imageData,
							width,
							height,
							width_step,
							_mask->widthStep,
							254,
							_est_bg_params ) )
	{
		//~ *err_stream << "Warning : Background pdf estimation!" << endl;
		return 1;
	}
	return 0;
}
	
int Mouloud :: estimate_radii ( 	double * radii,
									void * bg_params,
									void * obj_params,
									const double * data,
									const unsigned char * data_mask,
									unsigned int width,
									unsigned int height,
									unsigned int width_step,
									unsigned int mask_width_step )
{
	for ( unsigned int j = 0; j < width; ++ j )
	{
		double p_max = ((double*) _log_r_pdf->imageData)[j];
		radii[j] = 0;
		for ( unsigned int i = 1; i <= height; ++ i )
		{
			double v = ((double*) ( _log_r_pdf->imageData + i * _log_r_pdf->widthStep ) )[j % width];
			if ( v > p_max )
			{
				p_max = v;
				radii[j] = i;
			}
		}
	}
	
	
	
	
	return 0;
}
	
MOULOUD_COMPUTE_EDGES(unsigned char)
MOULOUD_COMPUTE_EDGES(char)
MOULOUD_COMPUTE_EDGES(unsigned short int)
MOULOUD_COMPUTE_EDGES(short int)
MOULOUD_COMPUTE_EDGES(long int)
MOULOUD_COMPUTE_EDGES(float)


