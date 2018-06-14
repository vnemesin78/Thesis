#include "c_iris_segmentation.hpp"

double r_ellipse ( const double & theta,
					const double * params )
{
	double a_2 = params[0] * params[0],
			b_2 = params[1] * params[1];
	return sqrt( 0.5 * ( a_2 + b_2 + ( a_2 - b_2 ) * cos ( 2 * ( theta + params[2] ) ) ) );
}


double d_ellipse ( const double & theta,
					const double * params )
{
	double a_2 = params[0] * params[0],
			b_2 = params[1] * params[1];
	return ( 0.5 * ( a_2 - b_2 ) * sin ( 2 * ( theta + params[2] ) ) / r_ellipse ( theta, params )  ) ;
}

double r_circle( 	const double & theta,
					const double * params )
{
	return params[0];
}


double d_circle ( const double & theta,
					const double * params )
{
	return 1 ;
}










c_iris_segmentation :: c_iris_segmentation ()
{
	initialize();
}

c_iris_segmentation :: c_iris_segmentation (	unsigned int width, //Largeur
												unsigned int height, //Hauteur
												unsigned int nb_directions,
												unsigned int kernel_size,
												unsigned int nb_samples,
												unsigned int dx,
												unsigned int nb_x,
												unsigned int dy,
												unsigned int nb_y,
												unsigned int r_min,
												unsigned int r_max,
												double r_ratio,
												unsigned int nb_r_ellipse,
												unsigned int dx_ellipse,
												unsigned int nb_x_ellipse,
												unsigned int dy_ellipse,
												unsigned int nb_y_ellipse,
												unsigned int nb_thetas,
												ostream * _err_stream
												 )
{
	initialize();
	if ( setup (	width,
					height,
					nb_directions,
					kernel_size,
					nb_samples,
					dx,
					nb_x,
					dy,
					nb_y,
					r_min,
					r_max,
					r_ratio,
					nb_r_ellipse,
					dx_ellipse,
					nb_x_ellipse,
					dy_ellipse,
					nb_y_ellipse,
					nb_thetas,
					_err_stream	) )
		throw ( invalid_argument("Arguments of c_iris_segmentation !") );

}

int c_iris_segmentation :: setup (	unsigned int width, //Largeur
									unsigned int height, //Hauteur
									unsigned int nb_directions,
									unsigned int kernel_size,
									unsigned int nb_samples,
									unsigned int dx,
									unsigned int nb_x,
									unsigned int dy,
									unsigned int nb_y,
									unsigned int r_min,
									unsigned int r_max,
									double r_ratio,
									unsigned int nb_r_ellipse,
									unsigned int dx_ellipse,
									unsigned int nb_x_ellipse,
									unsigned int dy_ellipse,
									unsigned int nb_y_ellipse,
									unsigned int nb_thetas,
									ostream * _err_stream
									 )
{
	free ();
	initialize();
	err_stream = _err_stream;
	if ( 	! width 		||
			! height 		||
			! nb_directions	||
			! nb_samples	||
			! nb_x			||
			! nb_y			||
			! dx			||
			! dy			||
			! r_min			||
			! r_max			||
			! r_ratio		||
			! nb_r_ellipse	||
			! nb_x_ellipse	||
			! nb_y_ellipse	||
			! nb_thetas		)
	{
		if ( err_stream )
			*err_stream << "Error: Invalid argument(s) in int c_iris_segmentation :: setup ( ... );" << endl;
		return 1;
	}

	_width = width;
	_height = height;
	_nb_points = nb_directions;
	_kernel_size = kernel_size;
	_nb_r = nb_samples;
	_dx = dx;
	_nb_x = nb_x;
	_dy = dy;
	_nb_y = nb_y;
	_r_min = r_min;
	_r_max = r_max;
	_r_ratio = r_ratio;
	_nb_r_ellipse = nb_r_ellipse;
	_dx_ellipse = dx_ellipse;
	_nb_x_ellipse = nb_x_ellipse;
	_dy_ellipse = dy_ellipse;
	_nb_y_ellipse = nb_y_ellipse;
	_nb_thetas = nb_thetas;
	if (integro_diff_ellipse_2.setup(	nb_directions,
										nb_samples ) )
		return 1;
	alloc();
	return 0;
}

int c_iris_segmentation :: setup( 	api_parameters & params,
										unsigned int width,
										unsigned int height,
										ostream * _err_stream,
										const char * prefix,
										const char * nb_directions_name,
										const char * kernel_size_name,
										const char * nb_samples_name,
										const char * dx_name,
										const char * nb_x_name,
										const char * dy_name,
										const char * nb_y_name,
										const char * r_min_name,
										const char * r_max_name,
										const char * r_ratio_name,
										const char * nb_r_ellipse_name,
										const char * dx_ellipse_name,
										const char * nb_x_ellipse_name,
										const char * dy_ellipse_name,
										const char * nb_y_ellipse_name,
										const char * nb_thetas_name
									)
{

	stringstream oss;
	int q = 0;
	err_stream = _err_stream;
	unsigned int nb_directions;
	{
		oss << prefix << "::" << nb_directions_name;
		if ( api_get_positive_integer( 	params,
										oss.str().c_str(),
										&nb_directions,
										err_stream ) )
			q = 1;
		oss.str("");
	}

	unsigned int kernel_size;
	{
		oss << prefix << "::" << kernel_size_name;
		if ( api_get_positive_integer( 	params,
										oss.str().c_str(),
										&kernel_size,
										err_stream ) )
			q = 1;
		oss.str("");
	}

	unsigned int nb_samples;
	{
		oss << prefix << "::" << nb_samples_name;
		if ( api_get_positive_integer( 	params,
										oss.str().c_str(),
										&nb_samples,
										err_stream ) )
			q = 1;
		oss.str("");
	}

	unsigned int dx;
	{
		oss << prefix << "::" << dx_name;
		if ( api_get_positive_integer( 	params,
										oss.str().c_str(),
										&dx,
										err_stream ) )
			q = 1;
		oss.str("");
	}

	unsigned int nb_x;
	{
		oss << prefix << "::" << nb_x_name;
		if ( api_get_positive_integer( 	params,
										oss.str().c_str(),
										&nb_x,
										err_stream ) )
			q = 1;
		oss.str("");
	}

	unsigned int dy;
	{
		oss << prefix << "::" << dy_name;
		if ( api_get_positive_integer( 	params,
										oss.str().c_str(),
										&dy,
										err_stream ) )
			q = 1;
		oss.str("");
	}

	unsigned int nb_y;
	{
		oss << prefix << "::" << nb_y_name;
		if ( api_get_positive_integer( 	params,
										oss.str().c_str(),
										&nb_y,
										err_stream ) )
			q = 1;
		oss.str("");
	}

	unsigned int r_min;
	{
		oss << prefix << "::" << r_min_name;
		if ( api_get_positive_integer( 	params,
										oss.str().c_str(),
										&r_min,
										err_stream ) )
			q = 1;
		oss.str("");
	}

	unsigned int r_max;
	{
		oss << prefix << "::" << r_max_name;
		if ( api_get_positive_integer( 	params,
										oss.str().c_str(),
										&r_max,
										err_stream ) )
			q = 1;
		oss.str("");
	}

	double r_ratio;
	{
		oss << prefix << "::" << r_ratio_name;
		if ( api_get_double( 	params,
								oss.str().c_str(),
								&r_ratio,
								err_stream ) )
			q = 1;
		oss.str("");
	}

	unsigned int nb_r_ellipse;
	{
		oss << prefix << "::" << nb_r_ellipse_name;
		if ( api_get_positive_integer( 	params,
										oss.str().c_str(),
										&nb_r_ellipse,
										err_stream ) )
			q = 1;
		oss.str("");
	}

	unsigned int dx_ellipse;
	{
		oss << prefix << "::" << dx_ellipse_name;
		if ( api_get_positive_integer( 	params,
										oss.str().c_str(),
										&dx_ellipse,
										err_stream ) )
			q = 1;
		oss.str("");
	}

	unsigned int nb_x_ellipse;
	{
		oss << prefix << "::" << nb_x_ellipse_name;
		if ( api_get_positive_integer( 	params,
										oss.str().c_str(),
										&nb_x_ellipse,
										err_stream ) )
			q = 1;
		oss.str("");
	}

	unsigned int dy_ellipse;
	{
		oss << prefix << "::" << dy_ellipse_name;
		if ( api_get_positive_integer( 	params,
										oss.str().c_str(),
										&dy_ellipse,
										err_stream ) )
			q = 1;
		oss.str("");
	}

	unsigned int nb_y_ellipse;
	{
		oss << prefix << "::" << nb_y_ellipse_name;
		if ( api_get_positive_integer( 	params,
										oss.str().c_str(),
										&nb_y_ellipse,
										err_stream ) )
			q = 1;
		oss.str("");
	}

	unsigned int nb_thetas;
	{
		oss << prefix << "::" << nb_thetas_name;
		if ( api_get_positive_integer( 	params,
										oss.str().c_str(),
										&nb_thetas,
										err_stream ) )
			q = 1;
		oss.str("");
	}

	if ( q )
		return -1;


	return (	setup (	width,
						height,
						nb_directions,
						kernel_size,
						nb_samples,
						dx,
						nb_x,
						dy,
						nb_y,
						r_min,
						r_max,
						r_ratio,
						nb_r_ellipse,
						dx_ellipse,
						nb_x_ellipse,
						dy_ellipse,
						nb_y_ellipse,
						nb_thetas,
						err_stream	)	);
}

int c_iris_segmentation :: segment_iris (	const IplImage * image,
												const IplImage * mask,
												const double & x_p,
												const double & y_p,
												double r_p )
{
	if (	! image 	||
			! mask 		)
	{
		if ( err_stream)
			*err_stream << "Error : Params in int c_iris_segmentation :: segment_iris (	... );" << endl;
		return -1;
	}
	if ( mask->depth != IPL_DEPTH_8U )
	{
		if ( err_stream )
			*err_stream << "Error : Mask format (not IPL_DEPTH_8U) in int c_iris_segmentation :: segment_iris (	... );" << endl;
		return -2;
	}

	//Détection d'un cercle
	{
		unsigned int nbs_values[3];
		polar_parameters min_v,
						   max_v,
						   opt_v;

		min_v.setup( 1 );
		max_v.setup( 1 );
		opt_v.setup( 1 );

		nbs_values[0] = _nb_x;
		nbs_values[1] = _nb_y;
		nbs_values[2] = _nb_r;

		min_v._x = x_p - _dx;
		min_v._y = y_p - _dy;
		min_v._params[0] = _r_min;
		if ( min_v._params[0] < r_p + _dx )
			min_v._params[0] = r_p + _dx;
		max_v._x = x_p + _dx;
		max_v._y = y_p + _dy;
		max_v._params[0] = _r_max;
		double v;
		if ( ( v = integro_diff_ellipse_2.compute( 	opt_v,
													image,
													mask,
													min_v,
													max_v,
													nbs_values,
													r_circle,
													d_circle,
													less_than,
													_kernel_size,
													_r_min,
													_r_max ) )	 == nan("") )
			return 1;
		if ( isnan(v ) )
			return 1;
		if ( v > -1.0e8 )
			return 1;
		_x_Daugman = opt_v._x;
		_y_Daugman = opt_v._y;
		_r_Daugman = opt_v._params[0];
	}

	//Debug

	//Affinage
	{
		unsigned int nbs_values[5];

		double radius_min,
				radius_max;

		polar_parameters min_v,
						   max_v,
						   opt_v;

		radius_min = (1 - _r_ratio) * _r_Daugman;
		if ( radius_min < _r_min )
			radius_min = _r_min;

		radius_max = (1 + _r_ratio) * _r_Daugman;
		if ( radius_max > _r_max )
			radius_max = _r_max;



		min_v.setup( 3 );
		max_v.setup( 3 );
		opt_v.setup( 3 );

		nbs_values[0] = _nb_x_ellipse;
		nbs_values[1] = _nb_y_ellipse;
		nbs_values[2] = _nb_r_ellipse;
		nbs_values[3] = _nb_r_ellipse;
		nbs_values[4] = _nb_thetas;

		min_v._x = _x_Daugman  - _dx_ellipse;
		min_v._y = _y_Daugman  - _dy_ellipse;
		min_v._params[0] = radius_min;
		min_v._params[1] = radius_min;
		min_v._params[2] = 0;

		max_v._x = _x_Daugman  + _dx_ellipse;
		max_v._y = _y_Daugman  + _dy_ellipse;
		max_v._params[0] = radius_max;
		max_v._params[1] = radius_max;
		max_v._params[2] = M_PI / 2 * ( _nb_thetas - 1.0 ) / _nb_thetas;
		double v;
		if ( ( v = integro_diff_ellipse_2.compute( 	opt_v,
													image,
													mask,
													min_v,
													max_v,
													nbs_values,
													r_ellipse,
													d_ellipse,
													less_than,
													_kernel_size,
													_r_min,
													_r_max ) )	 == nan("") )
			return 1;
		if ( isnan(v ) )
			return 1;

		_x_ellipse = opt_v._x;
		_y_ellipse = opt_v._y;
		_a_ellipse = opt_v._params[0];
		_b_ellipse = opt_v._params[1];
		_theta_ellipse = opt_v._params[2];
	}
	return 0;
}

int c_iris_segmentation :: segment_iris (	const IplImage * image,
												const IplImage * preprocessed_image,
												const IplImage * mask,
												const IplImage * s_mask,
												const double & x_p,
												const double & y_p,
												const double & a_p,
												const double & b_p,
												const double & theta_p )
{



	IplImage * _tmp_mask = cvCreateImage(	 cvSize( GET_IMAGE_WIDTH( preprocessed_image), GET_IMAGE_HEIGHT(preprocessed_image) ),
											IPL_DEPTH_8U,
											1 );

	memset( _tmp_mask->imageData, 255, _tmp_mask->widthStep * _tmp_mask->height );
	if ( segment_iris( 	preprocessed_image,
						_tmp_mask,
						x_p,
						y_p,
						MAX(a_p, b_p) ) )
	{
		cvReleaseImage( &_tmp_mask );
		return 1;
	}

	double f,
			dx,
			dy,
			r;

	if ( _b_ellipse > _a_ellipse )
	{
		r = _b_ellipse;
		f = _b_ellipse / _a_ellipse;
		dx = 1 + (f - 1) * cos( _theta_ellipse );
		dy = 1 + (f - 1) * sin( _theta_ellipse );
	}
	else
	{
		r = _a_ellipse;
		f = _a_ellipse / _b_ellipse;
		dx = 1 + (f - 1) * sin( _theta_ellipse );
		dy =  1 + (f - 1) * cos( _theta_ellipse );
	}
	cvSetImageROI ( _new_image,
					cvRect(	0, 0, GET_IMAGE_WIDTH( image ) * dx, GET_IMAGE_HEIGHT( image ) * dy ) );

	cvSetImageROI ( _new_mask,
					cvRect(	0, 0, GET_IMAGE_WIDTH( _new_image ), GET_IMAGE_HEIGHT( _new_image ) ) );

    dx = ( (double) GET_IMAGE_WIDTH( _new_image )) / GET_IMAGE_WIDTH( image);
    dy = ( (double) GET_IMAGE_HEIGHT( _new_image )) / GET_IMAGE_HEIGHT( image);



	image_resize ( 	(unsigned char*) _new_image->imageData,
					(unsigned char*) _new_mask->imageData,
					_new_image->roi->width,
					_new_image->roi->height,
					_new_image->widthStep,
					_new_mask->widthStep,
					(unsigned char*) get_image_data (image ),
					(unsigned char*) get_image_data (s_mask),
					GET_IMAGE_WIDTH( image ),
					GET_IMAGE_HEIGHT( image ),
					image->widthStep,
					s_mask->widthStep );
	cvReleaseImage( &_tmp_mask );


	_new_x = _x_ellipse * dx;
	_new_y = _y_ellipse * dy;
	_new_r = r;


	if ( dx == 1 && dy == 1 )
	{
		_new_x_pupil = x_p;
		_new_y_pupil = y_p;
		_new_a_pupil = a_p;
		_new_b_pupil = b_p;
		_new_theta_pupil = theta_p;
		
	}
	else
	{
		//Calcul des points de la pupille
		for ( unsigned int i = 0; i < _nb_points; ++ i )
		{
			double angle = 2 * M_PI / _nb_points * i;
			pupil_points[ 2 * i ] = ( x_p + cos ( theta_p ) * a_p * cos( angle ) + b_p * sin ( theta_p) * sin( angle )  );
			pupil_points[ 2 * i + 1 ] = ( y_p - sin ( theta_p ) * a_p * cos( angle ) + b_p * cos ( theta_p) * sin( angle )  );
			pupil_points[ 2 * i ] *= dx;
			pupil_points[ 2 * i + 1 ] *= dy;
		}
		
		unsigned char * _t_mask = new unsigned char[_nb_points];
		memset(_t_mask, 255, sizeof( char) * _nb_points );
		
		
		ellipse.fit( pupil_points, _nb_points, 1, _t_mask );
		delete[] _t_mask;
		
		_new_x_pupil = ellipse.x_center();
		_new_y_pupil = ellipse.y_center();
		_new_a_pupil = ellipse.a();
		_new_b_pupil = ellipse.b();
		_new_theta_pupil = ellipse.theta();
	}
	return 0;
}







c_iris_segmentation :: ~c_iris_segmentation()
{
	free();
	initialize();
}

void c_iris_segmentation :: alloc()
{
	_new_image = cvCreateImage(	cvSize( _width * (1 + 2 * _r_ratio), _height * ( 1 + 2 * _r_ratio)  ),
								IPL_DEPTH_8U,
								1 );
	_new_mask = cvCreateImage(	cvSize( _width * (1 + 2 * _r_ratio), _height * ( 1 + 2 * _r_ratio)  ),
								IPL_DEPTH_8U,
								1 );
	pupil_points = new double[ 2 * _nb_points];
	memset( pupil_points, 0, sizeof(double) * 2 * _nb_points );
}

void c_iris_segmentation :: free()
{
	if ( _new_image )
		cvReleaseImage (&_new_image );
	if ( _new_mask )
		cvReleaseImage (&_new_mask );
	if ( pupil_points)
		delete[] pupil_points;
}

void c_iris_segmentation :: initialize()
{
	_width = 0;
	_height = 0;
	_nb_points = 0;
	_kernel_size = 0;
	_dx = 0;
	_nb_x = 0;
	_dy = 0;
	_nb_y = 0;
	_r_min = 0;
	_r_max = 0;
	_nb_r = 0;
	_r_ratio = 0;
	_nb_r_ellipse = 0;
	_dx_ellipse = 0;
	_nb_x_ellipse = 0;
	_dy_ellipse = 0;
	_nb_y_ellipse = 0;
	_nb_thetas = 0;
	_x_Daugman = 0;
	_y_Daugman = 0;
	_r_Daugman = 0;
	 _x_ellipse = 0;
	_y_ellipse = 0;
	_a_ellipse = 0;
	_b_ellipse = 0;
	_theta_ellipse = 0;
	err_stream = NULL;
	_new_image = 0;
	_new_mask = 0;

	_new_x = 0;
	_new_y = 0;
	_new_r = 0;
	_new_x_pupil = 0;
	_new_y_pupil = 0;
	_new_a_pupil = 0;
	_new_b_pupil = 0;
	_new_theta_pupil = 0;
	pupil_points = 0;
}

