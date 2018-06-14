#include "c_integro_differential_2d.hpp"
#include "c_polar.hpp"
#include "filter.hpp"
#include <cmath>
#include <cstdlib> 
#include <complex> 
#include "image_utility.hpp"
#include <iostream>
using namespace std;
function_2d :: function_2d ( 	curve_function _x_curve,
								curve_function _y_curve )
{
	x_curve = _x_curve;
	y_curve = _y_curve;
}

void function_2d :: setup ( curve_function _x_curve,
							curve_function _y_curve)
{
	x_curve = _x_curve;
	y_curve = _y_curve;
}
void function_2d :: get_value (	double & x,
								double & y,
								const double & t,
								const double * params ) const
{
	x = x_curve( t, params );
	y = y_curve( t, params );
}

bool function_2d :: operator!() const
{
	return ( ! x_curve || ! y_curve );
}


c_integro_differential_2d :: c_integro_differential_2d ( )
{
	initialize();
}

c_integro_differential_2d :: c_integro_differential_2d (	unsigned int width,
															unsigned int height,
															unsigned int kernel_size )
{
	initialize();
	setup (	width,
			height,
			kernel_size );
}


int c_integro_differential_2d :: setup (	unsigned int width,
											unsigned int height,
											unsigned int kernel_size )
{
	free();
	initialize();
	
	//Checks
	if (	!width			||
			!height			||
			!kernel_size	)
	{
		*err_stream << "Error: Invalid argument(s) in int c_integro_differential_2d :: setup" << endl;
		return 1;
	}
	//INterpol
	{
		float f[4];
		f[0] = 1;
		f[1] = 1;
		f[2] = 1;
		f[3] = 1;
		interpol.setup(	width, 
						height,
						f );
	}
	
	img_32f = cvCreateImage(	cvSize(	width,
										height),
								IPL_DEPTH_32F,
								1 );
	h_grad = cvCreateImage(	cvSize(	width,
									height),
							IPL_DEPTH_32F,
							1 );
	v_grad = cvCreateImage(	cvSize(	width,
									height),
							IPL_DEPTH_32F,
							1 );
	
	_width = width;
	_height = height;
	_kernel_size = kernel_size;
	return 0;
}

double c_integro_differential_2d :: compute(	double * params,
												const IplImage * image,
												const IplImage * mask,
												const double * min_values,
												const double * max_values,
												unsigned int * nbs_values,
												unsigned int nb_params,
												const double & t_min,
												const double & t_max,
												unsigned int nb_t,
												const function_2d & f,
												const function_2d & df,
												order_function order )
{
	unsigned int 	image_width,
					image_height,
					mask_width,
					mask_height;
	
	//Check des arguments
	if ( 	!params			||
			!image			||
			!mask			||
			!min_values		||
			!max_values		||
			!nbs_values		||
			!nb_params		||
			!f 				||
			!df				||
			!order 			)
	{
		*err_stream << "Error : Invalid argument(s) in double c_integro_differential_par :: compute" << endl;
		return nan("");
	}
	for ( unsigned int i = 0; i < nb_params; ++ i )
	{
		if ( nbs_values[i] == 0 )
		{
			*err_stream << "Error : Invalid argument(s) in double c_integro_differential_par :: compute" << endl;
			return nan("");
		}
	}	
	image_width = GET_IMAGE_WIDTH(image);
	image_height = GET_IMAGE_HEIGHT(image);
	mask_width= GET_IMAGE_WIDTH(mask);
	mask_height = GET_IMAGE_HEIGHT(mask);
	
	if ( 	mask_width 	!= image_width	||
			mask_height != image_height	)
	{
		*err_stream << "Error : Different mask and image sizes in double c_integro_differential_par :: compute" << endl;
		return nan("");
	}
	if (	mask_width > _width		||
			mask_height > _height	)
	{
		*err_stream << "Error : Too large image sizes in double c_integro_differential_par :: compute" << endl;
		return nan("");
	}
	
	//ROIs
	cvSetImageROI( 	img_32f,
					cvRect(0, 0, image_width, image_height ) );
	cvSetImageROI( 	h_grad,
					cvRect(0, 0, image_width, image_height ) );
	cvSetImageROI( 	v_grad,
					cvRect(0, 0, image_width, image_height ) );
	
	//Interpolation
	if (	interpol.interpolate( 	image,
									mask,
									_kernel_size) )
	{
		*err_stream << "Error in double c_integro_differential_par :: compute" << endl;
		return nan("");
	}
	//Récupération de l'image
	interpol.get_image( img_32f );
	
	//Gradient
	cvSobel( 	img_32f, 
				h_grad, 
				1, 
				0, 
				_kernel_size );
				
	cvSobel( 	img_32f, 
				v_grad, 
				0, 
				1, 
				_kernel_size );
	
	//Maximisation de la coube
	double * curr_params = new double[nb_params];
	double max = search_maxima( params,
								curr_params,
								min_values,
								max_values,
								nbs_values,
								0,
								nb_params,
								mask,
								t_min,
								t_max,
								nb_t,
								f, 
								df,
								order );
	delete[] curr_params;
	return max;
}

double c_integro_differential_2d :: search_maxima( 	double * params,
													double * curr_params,
													const double * min_values,
													const double * max_values,
													unsigned int * nbs_values,
													unsigned int i,
													unsigned int nb_params,
													const IplImage * mask,
													const double & t_min,
													const double & t_max,
													unsigned int nb_t,
													const function_2d & f,
													const function_2d & df,
													order_function order )
{

	//Si i == nb_params, paramètres optimaux
	if ( i == nb_params )
	{
		memcpy (	params,
					curr_params,
					sizeof(double) * nb_params );
		return compute_integral( 	mask, 
									curr_params, 
									t_min,
									t_max,
									nb_t,
									f,
									df);
	}
	//Pas 
	double step,
		   * tmp_params = new double[nb_params],
		   max,
		   v;
	if ( nbs_values[i] == 1 )
		step = 0;
	else
		step = ( max_values[i] - min_values[i] ) / ( nbs_values[i] - 1 );
	//Max 0
	curr_params[i] = min_values[i];
	max = search_maxima( 	tmp_params,
							curr_params,
							min_values,
							max_values,
							nbs_values,
							i + 1,
							nb_params,
							mask,
							t_min,
							t_max,
							nb_t,
							f,
							df,
							order );
	memcpy (	params,
				tmp_params,
				sizeof(double) * nb_params );
				
	for ( unsigned int j = 1; j < nbs_values[i]; ++ j )
	{
		curr_params[i] = min_values[i] + j * step;
		v = search_maxima( 	tmp_params,
							curr_params,
							min_values,
							max_values,
							nbs_values,
							i + 1,
							nb_params,
							mask,
							t_min,
							t_max,
							nb_t,
							f,
							df,
							order );
		if ( order( v, max ) )
		{
			max = v;
			memcpy (	params,
						tmp_params,
						sizeof(double) * nb_params );
		}
	}
	
	delete[] tmp_params;
	return max;
	
}

double c_integro_differential_2d :: compute_integral( 	const IplImage * mask, 
														const double * curr_params, 
														const double & t_min,
														const double & t_max,
														unsigned int nb_t,
														const function_2d & f,
														const function_2d & df)
{
	//Données du masque
	unsigned int 	width,
					height,
					width_step,
					x_offset,
					y_offset;
	
	GET_IMAGE_DIM( 	mask, 
					width, 
					height, 
					width_step, 
					x_offset, 
					y_offset )
	
	
	double sum = 0;

	
	double t_step;
	if ( nb_t > 1 )
		t_step = ( t_max - t_min ) / (nb_t - 1);
	else
		t_step = 0;
	
	for ( double t = t_min; t < t_max; t += t_step )
	{
		double 	x_n,
				y_n,
				y_p,
				x_p;
		
		f.get_value( x_p, y_p, t, curr_params );
		df.get_value( x_n, y_n, t, curr_params );
		
		double v_h_grad,
			   v_v_grad;
			   
		unsigned char h_mask;
		interpol_2d ( 	v_h_grad,
						h_mask,
						x_p,
						y_p,
						(float*) h_grad->imageData,
						( (unsigned char*) ( mask->imageData + y_offset * width_step) ) + x_offset,
						width,
						height,
						h_grad->widthStep / sizeof(float),
						width_step );
						
		interpol_2d ( 	v_v_grad,
						h_mask,
						x_p,
						y_p,
						(float*) v_grad->imageData,
						( (unsigned char*) ( mask->imageData + y_offset * width_step) ) + x_offset,
						width,
						height,
						v_grad->widthStep / sizeof(float),
						width_step );
						
		//V normal
		if ( h_mask )
		{
			sum += x_n * v_v_grad; 
			sum += y_n * v_h_grad;
		}
	}
	return sum;
}

void c_integro_differential_2d :: initialize()
{
	_width = 0;
	_height = 0;
	_kernel_size = 0;
	err_stream = &cout;
	img_32f = NULL;
	h_grad = NULL;
	v_grad = NULL;
}

void c_integro_differential_2d :: free()
{
	if ( img_32f )
		cvReleaseImage ( &img_32f );
	if ( h_grad )
		cvReleaseImage ( &h_grad );
	if ( v_grad )
		cvReleaseImage ( &v_grad );
}

c_integro_differential_2d :: ~c_integro_differential_2d()
{
	free();
	initialize();
}
