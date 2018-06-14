#include "c_radon.hpp"
c_radon :: c_radon ( )
{
	initialize();
}
c_radon :: c_radon ( unsigned int r_size_max,
					 unsigned int nb_thetas_max )
{
	initialize();
	if ( setup ( 	r_size_max,
					nb_thetas_max ) )
		throw ( invalid_argument("Invalid arguments in c_radon :: c_radon ( unsigned int r_size_max,  unsigned int nb_thetas_max ); !") );
}

int c_radon :: setup (	unsigned int r_size_max,
						unsigned int nb_theta_max )
{
	free();
	initialize();
	
	if ( !r_size_max || !nb_theta_max )
	{
		*err_stream << "Error: Bad arguments in int c_radon :: setup (	unsigned int r_size_max, unsigned int nb_theta_max );" << endl;
		return 1;
	}
	_r_size_max = r_size_max;
	_nb_thetas_max = nb_theta_max;
	
	
	//~ _radon_transform = new double [ 2 * _r_size_max * _nb_thetas_max ];
	radon_image = cvCreateImage ( cvSize( _nb_thetas_max, 2 * _r_size_max  + 1),
								  IPL_DEPTH_64F,
								  1 );
	_radon_transform = (double*) radon_image->imageData;
	_width_step = radon_image->widthStep / sizeof(double);
	
	x_table = new double [ 2 * _r_size_max + 1];
	y_table = new double [ 2 * _r_size_max + 1];
	x_cos_table = new double [ 2 * _r_size_max + 1];
	y_sin_table = new double [ 2 * _r_size_max + 1];
	
	return 0;
}

template<class type> int c_radon :: compute_transform( 	const type * img_data, 
														unsigned int img_width,
														unsigned int img_height,
														unsigned int img_width_step,
														unsigned int nb_thetas )
{
	unsigned int r_size;
	unsigned int x_origin,
				 y_origin,
				 x_size,
				 y_size;
	double angle_step;
	if ( ! img_data || ! img_width || ! img_height || ! img_width_step || ! nb_thetas )
	{
		*err_stream << "Error: bad parameters in template<class type> int c_radon :: compute_transform( const type * img_data, unsigned int img_width, unsigned int img_height, unsigned int img_width_step, const double * thetas, unsigned int nb_thetas );" << endl;
		return -1;
	}
	r_size = (unsigned int) sqrt( img_width * img_width + img_height * img_height ) + 1;
	if ( r_size > _r_size_max )
	{
		*err_stream << "Error : too large image width or height!" << endl;
		return 1;
	}
	
	x_origin = img_width / 2;
	y_origin = img_height / 2;
	x_size = 2 * img_width;
	y_size = 2 * img_height;
	_nb_thetas = nb_thetas;
	_r_size = r_size;
	
	//Tables computation

	y_table[ y_size - 1 ] =  - 0.25 - y_origin;

	for ( unsigned int i = y_size - 2; i != (unsigned int) -1; -- i )
	{
		y_table[i] = y_table[i + 1] + 0.5;
	}

	x_table[ 0 ] =  - 0.25 - x_origin;
	for ( unsigned int i = 1; i < x_size; ++ i )
	{
		x_table[i] = x_table[i - 1] + 0.5;
	}

	// 0

	for ( unsigned int i = 0; i <= 2 * _r_size; ++ i )
	{
		memset ( 	_radon_transform + i * _nb_thetas_max,
					0,
					sizeof(double) * _nb_thetas );
	}

	angle_step = M_PI / nb_thetas;
	for ( unsigned int i = 0; i < nb_thetas; ++ i )
	{
		double angle = i * angle_step,
			   sinus = sin( angle ),
			   cosinus = cos( angle );
			   
		for ( unsigned int j = 0; j < y_size; ++ j )
			y_sin_table[j] = y_table[j] * sinus;
		
		for ( unsigned int j = 0; j < x_size; ++ j )
			x_cos_table[j] = x_table[j] * cosinus + r_size;

		for ( unsigned int j = 0; j < y_size; ++ j )
		{
			for ( unsigned int k = 0; k < x_size; ++ k )
			{

				unsigned int n_pixel = img_width_step * ( j / 2 ) + ( k / 2 );

				double pixel = img_data[n_pixel] /255;
				if ( pixel != 0 )
				{

					double 	r_ind,
							pixel_low;
					int r_low;
					pixel *= 0.25;
					r_ind = x_cos_table[k] + y_sin_table[j];
					r_low = (int) r_ind;
					pixel_low = pixel * ( 1.0 - r_ind + r_low );
					
					_radon_transform[ r_low * _width_step  + i ] += pixel_low;
					_radon_transform[ ( r_low + 1 ) * _width_step + i ] += pixel - pixel_low;
				}
				
			}
			
		}
	}

	
	
	//Debug
	//~ double max = 0;
	//~ for ( unsigned int i = 0; i < 2 * _r_size; ++ i )
	//~ {
		//~ for ( unsigned int j = 0; j < _nb_thetas; ++ j )
		//~ {
			//~ unsigned int n_pixel = i * _width_step + j;
			//~ if ( _radon_transform[n_pixel] > max )
				//~ max = _radon_transform[n_pixel];
		//~ }
	//~ 
	//~ }
	//~ 
	//~ 
	//~ for ( unsigned int i = 0; i < 2 * _r_size; ++ i )
	//~ {
		//~ for ( unsigned int j = 0; j < _nb_thetas; ++ j )
		//~ {
			//~ unsigned int n_pixel = i * _width_step + j;
			//~ _radon_transform[n_pixel] /= max;
		//~ }
	//~ 
	//~ }
	//~ 
	//~ 
	//~ 
	//~ 
	//~ 
	//~ 
	//~ cvSetImageROI( radon_image, cvRect (0,0, _nb_thetas, 2 * _r_size ) );
	//~ cvShowImage ( "ddd", radon_image );
	return 0;
}

template int c_radon :: compute_transform( 	const unsigned char * img_data, 
											unsigned int img_width,
											unsigned int img_height,
											unsigned int img_width_step,
											unsigned int nb_thetas );
											

template int c_radon :: compute_transform( 	const char * img_data, 
											unsigned int img_width,
											unsigned int img_height,
											unsigned int img_width_step,
											unsigned int nb_thetas );
											
template int c_radon :: compute_transform( 	const unsigned short int * img_data, 
											unsigned int img_width,
											unsigned int img_height,
											unsigned int img_width_step,
											unsigned int nb_thetas );
											

template int c_radon :: compute_transform( 	const short int * img_data, 
											unsigned int img_width,
											unsigned int img_height,
											unsigned int img_width_step,
											unsigned int nb_thetas );
											
template int c_radon :: compute_transform( 	const unsigned long int * img_data, 
											unsigned int img_width,
											unsigned int img_height,
											unsigned int img_width_step,
											unsigned int nb_thetas );
											

template int c_radon :: compute_transform( 	const long int * img_data, 
											unsigned int img_width,
											unsigned int img_height,
											unsigned int img_width_step,
											unsigned int nb_thetas );
											
template int c_radon :: compute_transform( 	const float * img_data, 
											unsigned int img_width,
											unsigned int img_height,
											unsigned int img_width_step,
											unsigned int nb_thetas );
											
template int c_radon :: compute_transform( 	const double * img_data, 
											unsigned int img_width,
											unsigned int img_height,
											unsigned int img_width_step,
											unsigned int nb_thetas );
											
											
int c_radon :: compute_transform ( 	const IplImage * image,
									unsigned int nb_thetas )
{
	if ( ! image )
	{
		*err_stream << "Error : Bad arguments in int c_radon :: compute_transform ( const IplImage * image, const double * thetas, unsigned int nb_theta );" << endl;
		return 1;
	}
	unsigned int width;
	unsigned int height;
	unsigned int width_step;
	
	switch( image->depth )
	{
		case ( IPL_DEPTH_8U ):
		{
			const unsigned char * image_data;
			width_step = image->widthStep;
			if ( image->roi != NULL )
			{
				image_data = ( (const unsigned char *) image->imageData ) + image->roi->yOffset * width_step + image->roi->xOffset;
				width = image->roi->width;
				height = image->roi->height;
			}
			else
			{
				image_data = ( (const unsigned char *) image->imageData );
				width = image->width;
				height = image->height;
			}
			return compute_transform(  image_data, width, height, width_step, nb_thetas );
		}
		break;
		case ( IPL_DEPTH_8S ):
		{
			const char * image_data;
			width_step = image->widthStep;
			if ( image->roi != NULL )
			{
				image_data = ( (const char *) image->imageData ) + image->roi->yOffset * width_step + image->roi->xOffset;
				width = image->roi->width;
				height = image->roi->height;
			}
			else
			{
				image_data = ( (const char *) image->imageData );
				width = image->width;
				height = image->height;
			}
			return compute_transform(  image_data, width, height, width_step, nb_thetas );
		}
		break;
		case ( IPL_DEPTH_16U ):
		{
			const unsigned short int * image_data;
			width_step = image->widthStep / sizeof( short int );
			if ( image->roi != NULL )
			{
				image_data = ( (const unsigned short int *) image->imageData ) + image->roi->yOffset * width_step + image->roi->xOffset;
				width = image->roi->width;
				height = image->roi->height;
			}
			else
			{
				image_data = ( (const unsigned short int *) image->imageData );
				width = image->width;
				height = image->height;
			}
			return compute_transform(  image_data, width, height, width_step, nb_thetas );
		}
		break;
		case ( IPL_DEPTH_16S ):
		{
			const short int * image_data;
			width_step = image->widthStep / sizeof( short int );
			if ( image->roi != NULL )
			{
				image_data = ( (const short int *) image->imageData ) + image->roi->yOffset * width_step + image->roi->xOffset;
				width = image->roi->width;
				height = image->roi->height;
			}
			else
			{
				image_data = ( (const short int *) image->imageData );
				width = image->width;
				height = image->height;
			}
			return compute_transform(  image_data, width, height, width_step, nb_thetas );
		}
		break;
		case ( IPL_DEPTH_32S ):
		{
			const long int * image_data;
			width_step = image->widthStep / sizeof( long int );
			if ( image->roi != NULL )
			{
				image_data = ( (const long int *) image->imageData ) + image->roi->yOffset * width_step + image->roi->xOffset;
				width = image->roi->width;
				height = image->roi->height;
			}
			else
			{
				image_data = ( (const long int *) image->imageData );
				width = image->width;
				height = image->height;
			}
			return compute_transform(  image_data, width, height, width_step, nb_thetas );
		}
		break;
		case ( IPL_DEPTH_32F ):
		{
			const float * image_data;
			width_step = image->widthStep / sizeof( float );
			if ( image->roi != NULL )
			{
				image_data = ( (const float *) image->imageData ) + image->roi->yOffset * width_step + image->roi->xOffset;
				width = image->roi->width;
				height = image->roi->height;
			}
			else
			{
				image_data = ( (const float *) image->imageData );
				width = image->width;
				height = image->height;
			}
			return compute_transform(  image_data, width, height, width_step, nb_thetas );
		}
		break;
		case ( IPL_DEPTH_64F ):
		{
			const double * image_data;
			width_step = image->widthStep / sizeof( double );
			if ( image->roi != NULL )
			{
				image_data = ( (const double *) image->imageData ) + image->roi->yOffset * width_step + image->roi->xOffset;
				width = image->roi->width;
				height = image->roi->height;
			}
			else
			{
				image_data = ( (const double *) image->imageData );
				width = image->width;
				height = image->height;
			}
			return compute_transform(  image_data, width, height, width_step, nb_thetas );
		}
		break;
		default:
			*err_stream << "Error : Unsupported image format in int c_radon :: compute_transform ( const IplImage * image, const double * thetas, unsigned int nb_theta );" << endl;
			return 1;
		break;
	}
	
}
											
c_radon :: ~c_radon()
{
	free();
	initialize();
}

void c_radon :: free()
{
	//~ if ( _radon_transform )
		//~ delete[] _radon_transform;
	if ( radon_image )
		cvReleaseImage ( &radon_image );
	if ( x_table )
		delete[] x_table;
	if ( y_table )
		delete[] y_table;
	if ( x_cos_table )
		delete[] x_cos_table;
	if ( y_sin_table )
		delete[] y_sin_table;
}

void c_radon :: initialize()
{
	_r_size_max = 0;
	_nb_thetas_max = 0;
	_radon_transform = 0;
	_r_size = 0;
	_nb_thetas = 0;
		
	err_stream = &cout;
	x_table = 0;
	y_table = 0;
	x_cos_table = 0;
	y_sin_table = 0;
	radon_image = 0;
}

double c_radon :: search_max ( 	double & r,
								double & theta ) const
{
	unsigned int size_r;
	int id_r = 0,
		id_theta = 0;
	double max;
	max = _radon_transform[0];
	size_r = _r_size * 2;
	for ( unsigned int i = 0; i < size_r; ++ i )
	{
		for ( unsigned int j = 0; j < _nb_thetas; ++ j )
		{
			double v = _radon_transform[ i * _width_step + j ];
			if ( v > max )
			{
				max = v;
				id_r = i;
				id_theta = j;
			}
			
		}
	}
	
	r = id_r - (double) _r_size + 0.5;
	theta = ( id_theta + 0.5 ) * M_PI / _nb_thetas; 
	return max;
}

