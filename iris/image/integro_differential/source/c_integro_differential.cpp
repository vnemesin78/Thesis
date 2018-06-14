#include "c_integro_differential.hpp"
#include "c_polar.hpp"
#include "filter.hpp"
#include <cmath>
#include <cstdlib> 
#include <complex> 
#include "image_utility.hpp"
#include <iostream>
using namespace std;
double compute_line_integral (	unsigned int & n_pixel,
								const double * polar_row,
								const unsigned char * mask_row,
								unsigned int first,
								unsigned int last )
{
	double value = 0;
	unsigned int i;
	for ( i = first; i < last; ++ i)
	{
		if ( mask_row[i] )
		{
			n_pixel ++;
			value += polar_row[i];
		}
	}
	return value;
}
c_integro_differential_operator :: c_integro_differential_operator()
{
	initialize();
}



c_integro_differential_operator :: c_integro_differential_operator ( 	unsigned int width,
																		unsigned int height,
																		unsigned int nb_directions,
																		unsigned int nb_samples,
																		unsigned int kernel_size )
{
	initialize();
	setup ( width,
			height,
			nb_directions,
			nb_samples,
			kernel_size );
}

void c_integro_differential_operator :: setup (	unsigned int width,
												unsigned int height,
												unsigned int nb_directions,
												unsigned int nb_samples,
												unsigned int kernel_size )
{
	unsigned int nb_dir;
	unsigned int dr_max;
	free();
	initialize();
	dr_max = nb_samples;
	nb_dir = nb_directions;
	polar_obj = new c_polar( 	dr_max, 
								nb_dir );

	_nb_directions = nb_directions;
	_nb_samples = nb_samples;
	
	polar_image = cvCreateImage ( 	cvSize( _nb_directions, _nb_samples ),
									IPL_DEPTH_64F,
									1 );
									
									
	_polar_data = ( double * ) polar_image->imageData;
	_polar_image_width_step = polar_image->widthStep / sizeof ( double );
	
	
	polar_mask = cvCreateImage ( 	cvSize( _nb_directions, _nb_samples ),
									IPL_DEPTH_8U,
									1 );
									
	_polar_mask = ( unsigned char * ) polar_mask->imageData;
	_polar_mask_width_step = polar_mask->widthStep;
	
	
	structuring_element_1 = cvCreateStructuringElementEx( 2 * kernel_size + 1,
														  2 * kernel_size + 1,
														  kernel_size,
														  kernel_size,
														  CV_SHAPE_ELLIPSE,
														  NULL);
	
	
	
	integrals = new double[nb_samples];
	memset( integrals,
			0,
			sizeof(double) * nb_samples );
			
	derivates = new double[nb_samples];
	memset( derivates,
			0,
			sizeof(double) * nb_samples );
			
	gaussian_kernel = new double[2 * kernel_size + 1];
	memset( gaussian_kernel,
			0,
			sizeof(double) * ( 2 * kernel_size + 1 ) );
			
	_kernel_size = kernel_size;
	
	//Noyau gaussian
	for (unsigned int n = 1; n <= kernel_size; ++ n)
	{
		gaussian_kernel[kernel_size + n] = exp( -(double)  (n * n) / 4 );
	}
	for (unsigned int n = 1; n <= kernel_size; ++ n)
	{
		gaussian_kernel[kernel_size - n] = exp( -(double)  (n * n) / 4 );
	}
	gaussian_kernel[ kernel_size ] = 1;
	
	
}


template <class type> double c_integro_differential_operator :: search_best_circle_params (	const type * data,
																								const unsigned char * mask,
																								unsigned int width,
																								unsigned int height,
																								unsigned int width_step,
																								unsigned int mask_width_step,
																								double x_min,
																								double x_max,
																								unsigned int n_x,
																								double y_min,
																								double y_max,
																								unsigned int n_y,
																								double r_min,
																								double r_max )
{

	
	
	double max = 0;
	_x = 0;
	_y = 0;
	_r = 0;
	
	//Test de validité
	if ( n_x < 2 || n_y < 2 )
		return nan("");
	
	for ( unsigned int i = 0; i < n_x; ++ i ) 
	{
		double x = x_min + i * (x_max - x_min ) / ( n_x - 1 );
		for ( unsigned int j = 0; j < n_y; ++ j ) 
		{
			double y = y_min + j * (y_max - y_min ) / ( n_y - 1 );
			unsigned int r = 0;
			//Transformée polaire
			polar_obj->compute( x, 
								y, 
								_nb_samples, 
								_nb_directions, 
								r_min, 
								r_max,  
								data,
								mask,
								width,
								height,
								width_step,
								mask_width_step );
								
			polar_obj->resize( 	_polar_data,
								_polar_mask,
								_nb_samples,
								_nb_directions,
								_polar_image_width_step,
								_polar_mask_width_step );
			//Suppression des zones erronnées (Mask)
			cvMorphologyEx(	polar_mask,
							polar_mask,
							NULL,
							structuring_element_1,
							CV_MOP_ERODE,
							1);
							
			compute_integral();
			compute_derivate();
			r = select_best_radius();
			if ( abs(derivates[r]) > max )
			{
				_x = x;
				_y = y;
				_r = ( (double) r ) / (_nb_samples - 1) * (r_max - r_min) + r_min;
				max = derivates[r];
			}

			
		}
	}
	return max;
	
	
}

double c_integro_differential_operator :: search_best_circle_params (	const IplImage * image,
																	const IplImage * mask,
																	double x_min,
																	double x_max,
																	unsigned int n_x,
																	double y_min,
																	double y_max,
																	unsigned int n_y,
																	double r_min,
																	double r_max )
{
	if ( image == 0 || mask == 0 )
		return nan("");
	if ( mask->depth != IPL_DEPTH_8U )
		return nan("");

	unsigned int width;
	unsigned int height;
	unsigned int width_step;
	//Mask
	const unsigned char * mask_data;
	unsigned int width_step_mask = mask->widthStep;
	if ( mask->roi )
		mask_data = (unsigned char*) mask->imageData + mask->roi->yOffset * width_step_mask + mask->roi->xOffset;
	else
		mask_data = (unsigned char*) mask->imageData;

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
			return search_best_circle_params (	image_data,
												mask_data,
												width,
												height,
												width_step,
												width_step_mask,
												x_min,
												x_max,
												n_x,
												y_min,
												y_max,
												n_y,
												r_min,
												r_max );
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
			return search_best_circle_params (	image_data,
												mask_data,
												width,
												height,
												width_step,
												width_step_mask,
												x_min,
												x_max,
												n_x,
												y_min,
												y_max,
												n_y,
												r_min,
												r_max );
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
			return search_best_circle_params (	image_data,
												mask_data,
												width,
												height,
												width_step,
												width_step_mask,
												x_min,
												x_max,
												n_x,
												y_min,
												y_max,
												n_y,
												r_min,
												r_max );
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
			return search_best_circle_params (	image_data,
												mask_data,
												width,
												height,
												width_step,
												width_step_mask,
												x_min,
												x_max,
												n_x,
												y_min,
												y_max,
												n_y,
												r_min,
												r_max );
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
			return search_best_circle_params (	image_data,
												mask_data,
												width,
												height,
												width_step,
												width_step_mask,
												x_min,
												x_max,
												n_x,
												y_min,
												y_max,
												n_y,
												r_min,
												r_max );
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
			return search_best_circle_params (	image_data,
												mask_data,
												width,
												height,
												width_step,
												width_step_mask,
												x_min,
												x_max,
												n_x,
												y_min,
												y_max,
												n_y,
												r_min,
												r_max );
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
			return search_best_circle_params (	image_data,
												mask_data,
												width,
												height,
												width_step,
												width_step_mask,
												x_min,
												x_max,
												n_x,
												y_min,
												y_max,
												n_y,
												r_min,
												r_max );
		}
		break;
		default:
			return nan("");
		break;
	}
}

c_integro_differential_operator :: ~c_integro_differential_operator()
{
	free();
	initialize();
}

void c_integro_differential_operator :: compute_integral()
{
	double * polar_row;
	unsigned char * polar_mask_row;
	for ( unsigned  i = 0; i < _nb_samples; ++ i)
	{
		unsigned int nb_pixels = 0;
		polar_row = _polar_data + i * _polar_image_width_step;
		polar_mask_row = _polar_mask + i * _polar_mask_width_step;
		
		//~ integrals[i] = compute_line_integral ( polar_row,
											   //~ polar_mask_row,
											   //~ 0,
											   //~ _nb_directions );
									   
		
		///Lignes suivantes inutiles vue la nature de l'opérateur intégro diff + détection "grossière"
		integrals[i] = compute_line_integral ( 	nb_pixels,
												polar_row,
												polar_mask_row ,
												0,
												_nb_directions / 8 );
											   
		integrals[i] += compute_line_integral (	nb_pixels,
												polar_row,
												polar_mask_row,
												( 3 * _nb_directions ) / 8,
												( 5 * _nb_directions ) / 8 );						   
										   
		integrals[i] += compute_line_integral (	nb_pixels,
												polar_row,
												polar_mask_row,
												( 7 * _nb_directions ) / 8,
												_nb_directions );
			integrals[i] /= 2 * ( ( double ) _nb_directions ) / ( nb_pixels  );
	}
	
}


void c_integro_differential_operator :: compute_derivate( )
{
	unsigned int end = _nb_samples - 1;
	for ( unsigned int i = 0; i < end; ++ i )
	{
		integrals[i] = integrals[ i + 1 ] - integrals[ i ];
	}
	
	//Application du noyeau gaussien
	filter_1d ( derivates,
				integrals,
				end,
				gaussian_kernel,
				2 * _kernel_size + 1 );
	

}

unsigned int c_integro_differential_operator :: select_best_radius( )
{
	double value = 0;
	unsigned int r;
	unsigned int end = _nb_samples - 1;
	//Recherche du max de la dérivée
	value = derivates[0];
	r = 0;
	for ( unsigned int i = 1; i < end; ++ i )
	{
		derivates[i] = abs( derivates[i] );
		if ( value < derivates[i] )
		{
			r = i;
		    value = derivates[i];
		}
	}
	return r;
}

void c_integro_differential_operator :: free()
{
	if ( polar_obj )
		delete polar_obj;
	if ( polar_image )
		cvReleaseImage ( &polar_image);
	if ( polar_mask )
		cvReleaseImage ( &polar_mask);
	if ( integrals )
		delete[] integrals;
	if ( derivates )
		delete[] derivates;
	if ( gaussian_kernel )
		delete[] gaussian_kernel;
	if ( structuring_element_1 )
		cvReleaseStructuringElement(&structuring_element_1);
}

void c_integro_differential_operator :: initialize()
{
	_x = 0; 
	_y = 0; 
	_r = 0;
					
			//Transformée polaire
	polar_obj = 0;
	_nb_directions = 0;
	_nb_samples = 0;
	_polar_data = 0;
	_polar_mask = 0;
	polar_mask = 0;
	polar_image = 0;
	structuring_element_1 = 0;
			//Noyau gaussian + maximisation
	integrals = 0;
	derivates = 0;
			
	gaussian_kernel = 0;
	_kernel_size = 0;
}

C_INTEGRODIFF_SEARCH_CIRCLE(unsigned char)
C_INTEGRODIFF_SEARCH_CIRCLE(char)
C_INTEGRODIFF_SEARCH_CIRCLE(unsigned short int)
C_INTEGRODIFF_SEARCH_CIRCLE(short int)
C_INTEGRODIFF_SEARCH_CIRCLE(unsigned long int)
C_INTEGRODIFF_SEARCH_CIRCLE(long int)
C_INTEGRODIFF_SEARCH_CIRCLE(float)
C_INTEGRODIFF_SEARCH_CIRCLE(double)
