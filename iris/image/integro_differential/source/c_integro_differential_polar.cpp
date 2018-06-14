#include "c_integro_differential_polar.hpp"
#include <iostream>
#include <cmath>
#include "image_utility.hpp"
using namespace std;
polar_parameters :: polar_parameters( void )
{
	_x = 0;
	_y = 0;
	_nb_params = 0;
	_params = NULL;
}

polar_parameters :: polar_parameters( unsigned int nb_params )
{
	_x = 0;
	_y = 0;
	setup ( nb_params );
}

polar_parameters & polar_parameters :: operator=( const polar_parameters & params )
{
	if ( _nb_params != params._nb_params )
	{
		setup ( params._nb_params );
	}
	_x = params._x;
	_y = params._y;
	memcpy( _params, 
			params._params,
			sizeof(double) * _nb_params );
	return *this;
}

ostream & operator << ( ostream & out, const polar_parameters & params )
{
	out << "x : " << params._x << endl;
	out << "y : " << params._y << endl;
	for ( unsigned int i = 0; i < params._nb_params; ++ i )
	{
		out << "params[" << i << "] : " << params._params[i] << endl; 
	}	
	return out;
}
						

void polar_parameters :: setup ( unsigned int nb_params )
{
	if ( _params )
		delete[] _params;
	_params = NULL;
		
	_nb_params = nb_params;
	if( _nb_params )
	{
		_params = new double [ _nb_params ];
		memset (	_params,
					0,
					sizeof( double ) * _nb_params );
	}
}

polar_parameters :: ~polar_parameters ()
{
	if ( _params )
		delete[] _params;
	_x = 0;
	_y = 0;
	_nb_params = 0;
	_params = NULL;
}


c_integro_differential_polar :: c_integro_differential_polar( void )
{
	initialize();
}

c_integro_differential_polar :: c_integro_differential_polar(	unsigned int nb_directions,
																		unsigned int nb_samples )
{
	initialize();
	setup ( nb_directions,
			nb_samples );
}

int c_integro_differential_polar :: setup ( 	unsigned int nb_directions,
												unsigned int nb_samples )
{
	free();
	initialize();
	
	if ( 	! nb_directions	||
			! nb_samples )
	{
		*err_stream << "Error in int c_integro_differential_polar :: setup : Invalid argument(s) !" << endl;
		return 1;
	}
	_nb_directions = nb_directions;
	_nb_samples = nb_samples;
	alloc();
	return 0;
}

c_integro_differential_polar :: ~c_integro_differential_polar( void )
{
	free();
	initialize();
}

double c_integro_differential_polar :: compute ( 	polar_parameters & computed_params,
														const IplImage * image,
														const IplImage * mask,
														const polar_parameters & min_values,
														const polar_parameters & max_values,
														unsigned int * nbs_values,
														curve_function radius_function,
														curve_function radius_derivate,
														order_function order,
														unsigned int kernel_size,
														const double & r_min,
														const double & r_max )
{
	if ( 	! image				||
			! mask				||
			! nbs_values		||
			! radius_function	||
			! radius_derivate	||
			! order				)
	{
		*err_stream << "Error in double c_integro_differential_polar :: compute : Invalid argument(s) !" << endl;
		*err_stream << 	"( 	! image				||"
						"! mask					||"
						"! nbs_values			||"
						"! radius_function		||"
						"! radius_derivate		||"
						"! order				)" << endl;
		return nan("");
	}
	
	if ( 	computed_params._nb_params != min_values._nb_params	||
			computed_params._nb_params != max_values._nb_params	)
	{
		*err_stream << "Error in double c_integro_differential_polar :: compute : Invalid argument(s) !" << endl;
		*err_stream << 	"( 	computed_params._nb_params != min_values._nb_params	||"
						"	computed_params._nb_params != max_values._nb_params	)" << endl;
		return nan("");
	}
	
	if ( 	GET_IMAGE_WIDTH( image ) != GET_IMAGE_WIDTH( mask )	||
			GET_IMAGE_HEIGHT( image ) != GET_IMAGE_HEIGHT( mask ) )
	{
		*err_stream << "Error in double c_integro_differential_polar :: compute : Invalid argument(s) !" << endl;
		*err_stream << 	"(	GET_IMAGE_WIDTH( image ) != GET_IMAGE_WIDTH( mask )	||"
						"	GET_IMAGE_HEIGHT( image ) != GET_IMAGE_HEIGHT( mask ) )" << endl;
		return nan("");	
	}
	
	if ( mask->depth != IPL_DEPTH_8U )
	{
		*err_stream << "Error in double c_integro_differential_polar :: compute : Invalid argument(s) !" << endl;
		*err_stream << 	"( mask->depth != IPL_DEPTH_8U )" << endl;
		return nan("");	
	}

	_nb_params = computed_params._nb_params;
	for ( unsigned int i = 0; i < _nb_params; ++ i )
	{
		if ( nbs_values[i] == 0 )
		{
			*err_stream << "Error in double c_integro_differential_polar :: compute : Invalid argument(s) !" << endl;
			*err_stream << 	"( _nb_values[" << i << "] == 0 )" << endl;
			return nan("");	
		}
	}
	
	_min_params = &min_values;
	_max_params = &max_values;
	_nbs_values = nbs_values;
	_radius_function = radius_function;
	_radius_derivate = radius_derivate;
	_order = order;
	_r_min = r_min;
	_r_max = r_max;
	
	_theta_step = 2 * M_PI / _nb_directions;
	_r_step = _nb_samples / (_r_max - _r_min );	
	//Pas
	_step.setup( _nb_params );
	if ( _nbs_values[0] > 1 )
		_step._x = ( _max_params->_x - _min_params->_x ) / ( _nbs_values[0] - 1 );
	if ( _nbs_values[1] > 1 )
		_step._y = ( _max_params->_y - _min_params->_y ) / ( _nbs_values[1] - 1 );	
	
	for ( unsigned int i = 0; i < _nb_params; ++ i )
	{
		if ( _nbs_values[i + 2] > 1 )
			_step._params[i] = ( _max_params->_params[i] - _min_params->_params[i]  ) / ( _nbs_values[i + 2] - 1 );	
	}
	//Alloc du tableau des paramètres
	if ( _nb_params )
		_params = new polar_parameters[_nb_params + 1];
	else
		_params = NULL;
	//Recherche du max
	double max = 0;

	computed_params = *_min_params;
	_params[0] = computed_params;	
	//Transformée polaire
	polar_obj->compute( _params[0]._x, 
						_params[0]._y, 
						_nb_samples, 
						_nb_directions, 
						r_min, 
						r_max,  
						image,
						mask );
	
	inter_2d->interpolate(	polar_obj->img(),
							polar_obj->p_mask() );
							
	inter_2d->get_image( img_32f );
	
	//Gradient
	cvSobel( 	img_32f, 
				_h_grad, 
				1, 
				0, 
				2 * kernel_size + 1 );
				
	cvSobel( 	img_32f, 
				_v_grad, 
				0, 
				1, 
				2 * kernel_size + 1 );

	max = compute_integral( computed_params );	
	for ( unsigned int x_id = 0; x_id < _nbs_values[0]; ++ x_id )
	{
		for ( unsigned int y_id = 0; y_id < _nbs_values[1]; ++ y_id )
		{

			//Transformée polaire
			polar_obj->compute( _params[0]._x, 
								_params[0]._y, 
								_nb_samples, 
								_nb_directions, 
								r_min, 
								r_max,  
								image,
								mask );
			
			inter_2d->interpolate(	polar_obj->img(),
									polar_obj->p_mask() );
									
			inter_2d->get_image( img_32f );
			
			//Gradient
			cvSobel( 	img_32f, 
						_h_grad, 
						1, 
						0, 
						2 * kernel_size + 1 );
						
			cvSobel( 	img_32f, 
						_v_grad, 
						0, 
						1, 
						2 * kernel_size + 1 );
			
			double v = search_maxima( 0 );
			if ( _order( v, max ) )
			{
				computed_params = _params[0];
				max = v;
			}
			
			_params[0]._y += _step._y;
		}
		_params[0]._x += _step._x;
		_params[0]._y = min_values._y;
	}
	if ( _params )
		delete[] _params;
	_params = NULL;
	
	return max;
}

void c_integro_differential_polar :: free()
{
	if ( polar_obj )
		delete polar_obj;
	if ( inter_2d )
		delete inter_2d;
	if ( _params )
		delete[] _params;
	if ( _h_grad )
		cvReleaseImage ( &_h_grad );
	if ( _v_grad )
		cvReleaseImage ( &_v_grad );
	if ( img_32f )
		cvReleaseImage ( & img_32f );
}

void c_integro_differential_polar :: initialize()
{
	_h_grad = 0;
	_v_grad = 0;
	img_32f = 0;
	err_stream = &cout;
	_nb_directions = 0;
	_nb_samples = 0;
	polar_obj = 0;
	inter_2d = 0;
	_min_params = 0;
	_max_params = 0;
	_nbs_values = 0;
	_nb_params = 0;
	_radius_derivate = 0;
	_radius_function = 0;
	_order = 0;
	_params = 0;
}

void c_integro_differential_polar :: alloc()
{
	_h_grad = cvCreateImage ( 	cvSize( _nb_directions, _nb_samples ),
								IPL_DEPTH_32F,
								1);
	
	_v_grad = cvCreateImage ( 	cvSize( _nb_directions, _nb_samples ),
								IPL_DEPTH_32F,
								1);
	
	img_32f = cvCreateImage ( 	cvSize( _nb_directions, _nb_samples ),
								IPL_DEPTH_32F,
								1);
	
	float factor[4];
		factor[0] = 1;
		factor[1] = 1;
		factor[2] = 1;
		factor[3] = 1;
		
	inter_2d = new c_interpol_2d ( 	_nb_directions, 
										_nb_samples,
										factor );
	
	polar_obj = new c_polar ( 	_nb_samples,
								_nb_directions );
	
}

double c_integro_differential_polar :: compute_integral ( const polar_parameters & params )
{
	double sum = 0;
	for ( unsigned int i = 0; i < _nb_directions; ++ i )
	{
		double theta, 
				x,
				r,
				y,
				dr,
				dx,
				v_v_grad = 0,
				v_h_grad = 0;
		unsigned char h_mask = 0;
		//Angle
		x = i;
		theta = i * _theta_step;	
		
		if ( theta < M_PI / 4 || ( theta > 3 * M_PI / 4 && theta < 5 * M_PI / 4) || ( theta > 7 * M_PI / 4 ) )
		{

			//Rayon de la courbe
			r = _radius_function( theta , 
								  params._params );
			y = ( r - _r_min ) * _r_step;
		
		
			//Derivée	
			dr = _radius_derivate(	theta, 
									params._params ); 
			dx = r / _r_step;
			

			
			
			
			interpol_2d ( 	v_h_grad,
							h_mask,
							x,
							y,
							(float*) _h_grad->imageData,
							polar_obj->image_mask(),
							_nb_directions,
							_nb_samples,
							_h_grad->widthStep / sizeof(float),
							polar_obj->mask_width_step() );
							
			interpol_2d ( 	v_v_grad,
							h_mask,
							x,
							y,
							(float*) _v_grad->imageData,
							polar_obj->image_mask(),
							_nb_directions,
							_nb_samples,
							_v_grad->widthStep / sizeof(float),
							polar_obj->mask_width_step() );
						
			if ( h_mask )
			{
				double norm = sqrt( dx *dx + dr * dr );
				
				sum += - dx * v_v_grad / norm; 
				sum += dr * v_h_grad / norm;
			}
		}
	}

	return sum;
}

double c_integro_differential_polar :: search_maxima( unsigned int i )
{
	double max;
	if ( i == _nb_params )
		return compute_integral( _params[i] );
		
	_params[i]._params[i] = _min_params->_params[i];
	_params[i + 1] = _params[i];
	
	max = search_maxima ( i + 1 );
	_params[i] = _params[i + 1];
	
	for ( unsigned int id = 1; id < _nbs_values[ i + 2 ]; ++ id )
	{
		_params[i + 1]._params[i] += _step._params[i];	
		double v = search_maxima ( i + 1 );
		if ( _order ( v, max ) )
		{
			max = v;
			_params[i] = _params[i + 1];			
		}
	}
	return max;
}





