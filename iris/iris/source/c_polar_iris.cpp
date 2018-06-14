#include "c_polar_iris.hpp"

c_polar_iris :: c_polar_iris( void )
{
	initialize();
}

c_polar_iris :: c_polar_iris(	unsigned int nb_directions,
									unsigned int nb_samples,
									unsigned int nb_samples_iris,
									unsigned int width,
									unsigned int height,
										ostream * _err_stream )
{
	initialize();
	if ( setup ( 	nb_directions,
					nb_samples,
					nb_samples_iris,
					width,
					height,
					_err_stream ) )
		throw ( invalid_argument ( "Argument(s) of c_polar_iris(...)") );
}
int c_polar_iris :: setup (	unsigned int nb_directions,
							unsigned int nb_samples,
							unsigned int nb_samples_iris,
							unsigned int width,
							unsigned int height,
							ostream * _err_stream )
{
	unsigned int nb_samples_0 = (unsigned int) sqrt( width * width + height * height ) + 1;
	free();
	initialize();
	err_stream = _err_stream;
	if ( ! nb_directions 	||
		 ! nb_samples 		||
		 ! width			||
		 ! height			||
		 ! nb_samples_iris	)
	{
		if ( err_stream )
			(*err_stream) << "Error : Argument(s) of c_polar_iris :: setup !" << endl;
		return 1;
	}

	_nb_directions = nb_directions;
	_nb_samples = nb_samples;
	_nb_samples_iris = nb_samples_iris;
	polar_obj = new c_polar( nb_samples_0,
							 (unsigned int) 4 * M_PI * ( 2 * sqrt( width * width + height * height ) + 2.5 ) );

	_tmp_polar_image = cvCreateImage ( 	cvSize ( 	nb_directions,
													nb_samples_0 ),
										IPL_DEPTH_64F,
										1 );

	_tmp_polar_mask = cvCreateImage ( 	cvSize ( 	nb_directions,
													nb_samples_0 ),
										IPL_DEPTH_8U,
										1 );

	_polar_image = cvCreateImage ( cvSize (	nb_directions,
											nb_samples ),
									IPL_DEPTH_64F,
									1 );
	_polar_mask = cvCreateImage ( cvSize (	nb_directions,
											nb_samples ),
									IPL_DEPTH_8U,
									1 );
	_pupil_radii = new double[_nb_directions];
	_f_radii = new double[_nb_directions];
	_iris_radii = new double[_nb_directions];

	return 0;
}

int c_polar_iris :: setup ( 	api_parameters & params,
								unsigned int width,
								unsigned int height,
								ostream * _err_stream,
								const char * n_space,
								const char * nb_dir_name,
								const char * nb_samples_name,
								const char * nb_samples_iris_name )
{
	int q = 0;
	unsigned int 	nb_directions,
					nb_samples,
					nb_samples_iris;

	stringstream oss;
	err_stream = _err_stream;
	oss << n_space << "::" << nb_dir_name;
	if ( api_get_positive_integer( params,
								   oss.str().c_str(),
								   &nb_directions,
								   err_stream ) )
		q = 1;
	oss.str("");

	oss << n_space << "::" << nb_samples_name;
	if ( api_get_positive_integer( params,
								   oss.str().c_str(),
								   &nb_samples,
								   err_stream ) )
		q = 1;
	oss.str("");

	oss << n_space << "::" << nb_samples_iris_name;
	if ( api_get_positive_integer( params,
								   oss.str().c_str(),
								   &nb_samples_iris,
								   err_stream ) )
		q = 1;
	oss.str("");
	if ( q )
		return 1;

	return ( setup( 	nb_directions,
						nb_samples,
						nb_samples_iris,
						width,
						height,
						err_stream ) );
}


int c_polar_iris :: compute (	const IplImage * image,
								const IplImage * mask,
								const double & x_p,
								const double & y_p,
								const double & r_p,
								const double & x_i,
								const double & y_i,
								const double & r_i )
{
	return compute ( 	image,
						mask,
						x_p,
						y_p,
						r_p,
						r_p,
						0,
						x_i,
						y_i,
						r_i );
}





int c_polar_iris :: compute (	const IplImage * image,
								const IplImage * mask,
								const double & x_p,
								const double & y_p,
								const double & a_p,
								const double & b_p,
								const double & theta_p,
								const double & x_i,
								const double & y_i,
								const double & r_i )
{
	if ( ! image || ! mask )
	{
		//*err_stream << "Error : No Image or mask found!" << endl;
		return -1;
	}

	unsigned int 	img_width,
					img_height,
					img_x_offset,
					img_y_offset,
					img_width_step,
					mask_width,
					mask_height,
					mask_x_offset,
					mask_y_offset,
					mask_width_step;

	GET_IMAGE_DIM( image,
				   img_width,
				   img_height,
				   img_width_step,
				   img_x_offset,
				   img_y_offset );

	GET_IMAGE_DIM( mask,
				   mask_width,
				   mask_height,
				   mask_width_step,
				   mask_x_offset,
				   mask_y_offset );

	if ( 	img_width != mask_width 		||
			img_height != mask_height 		||
			mask->depth != IPL_DEPTH_8U		)
	{
		if ( err_stream)
		{
			*err_stream << "Error : Image(s) dimension(s)!" << endl;
			*err_stream << "width" << endl;
			*err_stream << img_width << "\t" << mask_width << "\t"<< endl;
			*err_stream << "height" << endl;
			*err_stream << img_height << "\t" << mask_height << "\t"<< endl;
		}
		return 1;
	}

	switch (image->depth)
	{
		case( IPL_DEPTH_8S):
			return ( compute( 	(char * ) ( image->imageData + img_width_step * img_y_offset ) + img_x_offset,
								(unsigned char *) mask->imageData + mask_x_offset + mask_width_step * mask_y_offset,
								img_width,
								img_height,
								img_width_step,
								mask_width_step,
								x_p,
								y_p,
								a_p,
								b_p,
								theta_p,
								x_i,
								y_i,
								r_i ) );
		break;
		case( IPL_DEPTH_8U):
			return ( compute( 	(unsigned char * ) ( image->imageData + img_width_step * img_y_offset ) + img_x_offset,
								(unsigned char *) mask->imageData + mask_x_offset + mask_width_step * mask_y_offset,
								img_width,
								img_height,
								img_width_step,
								mask_width_step,
								x_p,
								y_p,
								a_p,
								b_p,
								theta_p,
								x_i,
								y_i,
								r_i ) );
		break;
		case( IPL_DEPTH_16S):
			return ( compute( 	(short int * ) ( image->imageData + img_width_step * img_y_offset ) + img_x_offset,
								(unsigned char *) mask->imageData + mask_x_offset + mask_width_step * mask_y_offset,
								img_width,
								img_height,
								img_width_step / sizeof( short int),
								mask_width_step,
								x_p,
								y_p,
								a_p,
								b_p,
								theta_p,
								x_i,
								y_i,
								r_i ) );
		break;
		case( IPL_DEPTH_16U):
			return ( compute( 	(unsigned short int * ) ( image->imageData + img_width_step * img_y_offset ) + img_x_offset,
								(unsigned char *) mask->imageData + mask_x_offset + mask_width_step * mask_y_offset,
								img_width,
								img_height,
								img_width_step / sizeof( short int),
								mask_width_step,
								x_p,
								y_p,
								a_p,
								b_p,
								theta_p,
								x_i,
								y_i,
								r_i ) );
		break;
		case( IPL_DEPTH_32S):
			return ( compute( 	(long int * ) ( image->imageData + img_width_step * img_y_offset ) + img_x_offset,
								(unsigned char *) mask->imageData + mask_x_offset + mask_width_step * mask_y_offset,
								img_width,
								img_height,
								img_width_step / sizeof(long int),
								mask_width_step,
								x_p,
								y_p,
								a_p,
								b_p,
								theta_p,
								x_i,
								y_i,
								r_i ) );
		break;
		case( IPL_DEPTH_32F):
			return ( compute( 	(float * ) ( image->imageData + img_width_step * img_y_offset ) + img_x_offset,
								(unsigned char *) mask->imageData + mask_x_offset + mask_width_step * mask_y_offset,
								img_width,
								img_height,
								img_width_step / sizeof(float),
								mask_width_step,
								x_p,
								y_p,
								a_p,
								b_p,
								theta_p,
								x_i,
								y_i,
								r_i ) );
		break;
		case( IPL_DEPTH_64F):
			return ( compute( 	(double * ) ( image->imageData + img_width_step * img_y_offset ) + img_x_offset,
								(unsigned char *) mask->imageData + mask_x_offset + mask_width_step * mask_y_offset,
								img_width,
								img_height,
								img_width_step / sizeof(double),
								mask_width_step,
								x_p,
								y_p,
								a_p,
								b_p,
								theta_p,
								x_i,
								y_i,
								r_i ) );
		break;
		default:
			if ( err_stream )
				*err_stream << "Error : Unsupported image format!" << endl;
			return -2;
		break;
	}
}

template <class type> int c_polar_iris :: compute( 	const type * img_data,
														const unsigned char * iris_mask,
														unsigned int width,
														unsigned int height,
														unsigned int img_width_step,
														unsigned int mask_width_step,
														const double & x_p,
														const double & y_p,
														const double & a_p,
														const double & b_p,
														const double & theta_p,
														const double & x_i,
														const double & y_i,
														const double & r_i )
{
	double r_min,
			r_max;

	unsigned int 	nb_samples_0,
					nb_directions_0;

	double dy,
			y,
			y_step,
			y_end;

	//Calcul des rayons
	compute_pupil_radii( x_p,
						 y_p,
						 a_p,
						 b_p,
						 theta_p );
	compute_iris_radii( x_i,
						y_i,
						r_i,
						x_p,
						y_p );
	compute_f_radii	( );


	//Recherche des rayons extrèmes
	r_min = _pupil_radii[0];
	r_max = _f_radii[0];
	for ( unsigned int i = 0; i < _nb_directions; ++ i )
	{

		if ( _pupil_radii[i] < r_min )
			r_min = _pupil_radii[i];
		if ( _f_radii[i] > r_max )
			r_max = _f_radii[i];
	}

	//Paramètres pour la transformée polaire
	nb_samples_0 = (unsigned int) 2 * ( r_max - r_min + 1 );
	nb_directions_0 = (unsigned int) 4 * M_PI * 2 * r_max + 1;
	y_step = nb_samples_0 / ( r_max - r_min );

	//Transformée polaire
	if ( polar_obj->compute ( 	x_p,
								y_p,
								nb_samples_0,
								nb_directions_0,
								r_min,
								r_max,
								img_data,
								iris_mask,
								width,
								height,
								img_width_step,
								mask_width_step ) )
	{
		if (err_stream )
			*err_stream << "Polar transform!" << endl;
		return 1;
	}
	//Redim.
	cvSetImageROI ( _tmp_polar_image,
					cvRect( 0,
							0,
							_nb_directions,
							nb_samples_0 ) );
	cvSetImageROI ( _tmp_polar_mask,
					cvRect( 0,
							0,
							_nb_directions,
							nb_samples_0 ) );
	if ( polar_obj->resize( _tmp_polar_image,
							_tmp_polar_mask ) )
	{

		if (err_stream )
			*err_stream << "Polar transform resizing!" << endl;
		return 1;
	}
	for ( unsigned int i = 0; i < _nb_directions; ++ i)
	{
		//Debut
		y = ( _pupil_radii[i] - r_min ) * y_step;

		//Fin
		y_end = ( _f_radii[i] - r_min ) * y_step;

		//Pas
		dy = ( y_end - y ) / _nb_samples;
		if ( dy > 0 )
		{
			for ( unsigned int j = 0; j < _nb_samples; ++ j )
			{
				resample_2d(	( (double *) ( _polar_image->imageData + j * _polar_image->widthStep ) )[i],
								( (unsigned char *) ( _polar_mask->imageData + j * _polar_mask->widthStep ) )[i],
								i,
								y,
								1,
								dy,
								(double *) ( _tmp_polar_image->imageData ),
								(unsigned char *) ( _tmp_polar_mask->imageData),
								_nb_directions,
								nb_samples_0,
								_polar_image->widthStep / sizeof(double),
								_polar_mask->widthStep  );

				y += dy;
				//~ ( (double *) ( _polar_image->imageData + j * _polar_image->widthStep ) )[i] /= 255.0;
			}

		}
		else
		{
			if ( err_stream )
				*err_stream << "Error : Collision between iris and pupil!" << endl;
			return 1;
		}
	}
	


	
	
	
	
	return 0;

}

c_polar_iris :: ~c_polar_iris()
{
	free();
	initialize();
}


void c_polar_iris :: free()
{
	if ( polar_obj )
		delete polar_obj;

	if ( _tmp_polar_image )
		cvReleaseImage ( &_tmp_polar_image );

	if ( _tmp_polar_mask )
		cvReleaseImage ( &_tmp_polar_mask );

	if ( _polar_image )
		cvReleaseImage ( &_polar_image );

	if ( _polar_mask )
		cvReleaseImage ( &_polar_mask );

	if ( _pupil_radii )
		delete[] _pupil_radii;
	if ( _iris_radii )
		delete[] _iris_radii;
	if ( _f_radii )
		delete[] _f_radii;
}


void c_polar_iris :: initialize()
{
	_nb_samples_iris = 0;
	_nb_directions = 0;
	_nb_samples = 0;
			//PT
	polar_obj = 0;

	_tmp_polar_image = 0;
	_tmp_polar_mask = 0;

	_polar_image = 0;
	_polar_mask = 0;

	err_stream = NULL;
	_pupil_radii = 0;
	_iris_radii	= 0;
	_f_radii = 0;
}

void c_polar_iris :: compute_pupil_radii( 	const double & x_p,
												const double & y_p,
												const double & a_p,
												const double & b_p,
												const double & theta_p )
{
	double theta_step = 2 * M_PI / _nb_directions,
			a_2 = a_p * a_p,
			b_2 = b_p * b_p;
	for ( unsigned int i = 0; i < _nb_directions; ++ i )
	{
		_pupil_radii[i] = sqrt( 0.5 * ( a_2 + b_2 + ( a_2 - b_2 ) * cos ( 2 * ( i * theta_step - theta_p ) ) ) );
	}
}

void c_polar_iris :: compute_iris_radii(	const double & x_i,
											const double & y_i,
											const double & r_i,
											const double & x_p,
											const double & y_p )
{
	double theta_step = 2 * M_PI / _nb_directions,
			dx = x_i - x_p,
			dy = y_i - y_p,
			d_2 = dx * dx + dy * dy,
			d = sqrt(d_2),
			phi = atan2( dy, dx ),
			r_2 = r_i * r_i;

	double theta,
			zeta,
			cos_zeta,
			delta;
	for ( unsigned int i = 0; i < _nb_directions; ++ i )
	{
		theta = i * theta_step;
		zeta = phi - theta;
		cos_zeta = cos( zeta );
		delta = d_2 * cos_zeta * cos_zeta + ( r_2 - d_2 );
		_iris_radii[i] = d * cos_zeta + sqrt(delta);
	}
}

void c_polar_iris :: compute_f_radii( )
{
	double f = ( _nb_samples - _nb_samples_iris )  / ( (double) _nb_samples_iris );
	for ( unsigned int i = 0; i < _nb_directions; ++ i )
	{
		_f_radii[i] = _iris_radii[i] + f * ( _iris_radii[i] - _pupil_radii[i] );
	}
}

















C_POLAR_IRIS_COMPUTE(char)
C_POLAR_IRIS_COMPUTE(unsigned char)
C_POLAR_IRIS_COMPUTE(short int)
C_POLAR_IRIS_COMPUTE(unsigned short int)
C_POLAR_IRIS_COMPUTE(long int)
C_POLAR_IRIS_COMPUTE(unsigned long int)
C_POLAR_IRIS_COMPUTE(float)
C_POLAR_IRIS_COMPUTE(double)
