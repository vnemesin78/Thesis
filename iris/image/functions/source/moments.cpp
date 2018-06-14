#include "moments.hpp"
#include <iostream>
using namespace std;
template<class type> double get_image_mean ( const type * image_data,
											 unsigned int width,
											 unsigned int height,
											 unsigned int width_step )
{
	double sum = 0; 
	for ( unsigned int i = 0; i < height; ++ i )
	{
		for ( unsigned int j = 0; j < width; ++ j )
		{
			sum += image_data[ i * width_step + j ];
		}
	}
	sum /= ( width * height );
	return sum;
}

GET_IMAGE_MEAN_TEMP( unsigned char )
GET_IMAGE_MEAN_TEMP( char )

GET_IMAGE_MEAN_TEMP( unsigned short int )
GET_IMAGE_MEAN_TEMP( short int )

GET_IMAGE_MEAN_TEMP( unsigned long int )
GET_IMAGE_MEAN_TEMP( long int )

GET_IMAGE_MEAN_TEMP( float )
GET_IMAGE_MEAN_TEMP( double )

double get_image_mean( const IplImage * image )
{
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
			return get_image_mean ( image_data, width, height, width_step );
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
			return get_image_mean ( image_data, width, height, width_step );
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
			return get_image_mean ( image_data, width, height, width_step );
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
			return get_image_mean ( image_data, width, height, width_step );
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
			return get_image_mean ( image_data, width, height, width_step );
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
			return get_image_mean ( image_data, width, height, width_step );
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
			return get_image_mean ( image_data, width, height, width_step );
		}
		break;
		default:
			return 0;
		break;
	}
}

template<class type> double get_image_variance ( 	const type * image_data,
													unsigned int width,
													unsigned int height,
													unsigned int width_step,
													double mean )
{
	double 	sum = 0,
			tmp;
	for ( unsigned int i = 0; i < height; ++ i )
	{
		for ( unsigned int j = 0; j < width; ++ j )
		{
			tmp = image_data[ i * width_step + j ] - mean;
			sum += tmp * tmp;
		}
	}
	sum /= ( width * height );
	return sum;
}

GET_IMAGE_VARIANCE_TEMP( unsigned char )
GET_IMAGE_VARIANCE_TEMP( char )

GET_IMAGE_VARIANCE_TEMP( unsigned short int )
GET_IMAGE_VARIANCE_TEMP( short int )

GET_IMAGE_VARIANCE_TEMP( unsigned long int )
GET_IMAGE_VARIANCE_TEMP( long int )

GET_IMAGE_VARIANCE_TEMP( float )
GET_IMAGE_VARIANCE_TEMP( double )


double get_image_variance( 	const IplImage * image,
							double mean  )
{
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
			return get_image_variance( image_data, width, height, width_step, mean );
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
			return get_image_variance( image_data, width, height, width_step, mean );
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
			return get_image_variance( image_data, width, height, width_step, mean );
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
			return get_image_variance( image_data, width, height, width_step, mean );
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
			return get_image_variance( image_data, width, height, width_step, mean );
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
			return get_image_variance( image_data, width, height, width_step, mean );
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
			return get_image_variance( image_data, width, height, width_step, mean );
		}
		break;
		default:
			return 0;
		break;
	}
}

template<class type> double get_image_mean ( const type * image_data,
											 const unsigned char * mask_data,
											 unsigned int width,
											 unsigned int height,
											 unsigned int width_step,
											 unsigned int mask_width_step,
											 unsigned char mask_value )
{
	double sum = 0;
	unsigned int nb_pixels = 0; 
	for ( unsigned int i = 0; i < height; ++ i )
	{
		for ( unsigned int j = 0; j < width; ++ j )
		{
			if ( mask_data[ i * mask_width_step + j] == mask_value )
			{
				sum += image_data[ i * width_step + j ];
				++ nb_pixels;
			}
		}
	}
	if ( nb_pixels == 0 )
		return nan("");
	sum /= ( nb_pixels );
	return sum;
}

GET_IMAGE_MEAN_TEMP_2( unsigned char )
GET_IMAGE_MEAN_TEMP_2( char )

GET_IMAGE_MEAN_TEMP_2( unsigned short int )
GET_IMAGE_MEAN_TEMP_2( short int )

GET_IMAGE_MEAN_TEMP_2( unsigned long int )
GET_IMAGE_MEAN_TEMP_2( long int )

GET_IMAGE_MEAN_TEMP_2( float )
GET_IMAGE_MEAN_TEMP_2( double )

double get_image_mean( 	const IplImage * image,
						const IplImage * mask,
						unsigned char mask_value  )
{
	unsigned int width,
				 width_bis,
				 height_bis,
				 height;
	unsigned int width_step;
	if ( mask->depth != IPL_DEPTH_8U )
		return 0;
	
	
	if ( image->roi != NULL )
	{
		width = image->roi->width;
		height = image->roi->height;
	}
	else
	{
		width = image->width;
		height = image->height;
	}
	
	const unsigned char * mask_data;
	unsigned int mask_width_step = mask->widthStep;
	if ( mask->roi != NULL )
	{
		width_bis = mask->roi->width;
		height_bis = mask->roi->height;
		mask_data = ( const unsigned char * ) mask->imageData + mask->roi->yOffset * mask_width_step + mask->roi->xOffset;
	}
	else
	{
		width_bis = mask->width;
		height_bis = mask->height;
		mask_data = ( const unsigned char * ) mask->imageData ;
	}
	
	if ( width_bis != width || height_bis != height )
		return 0;
	
	switch( image->depth )
	{
		case ( IPL_DEPTH_8U ):
		{
			const unsigned char * image_data;
			width_step = image->widthStep;
			if ( image->roi != NULL )
			{
				image_data = ( (const unsigned char *) image->imageData ) + image->roi->yOffset * width_step + image->roi->xOffset;

			}
			else
			{
				image_data = ( (const unsigned char *) image->imageData );
			}
			return get_image_mean ( image_data,
									mask_data,
									width, 
									height, 
									width_step,
									mask_width_step,
									mask_value );
		}
		break;
		case ( IPL_DEPTH_8S ):
		{
			const char * image_data;
			width_step = image->widthStep;
			if ( image->roi != NULL )
			{
				image_data = ( (const char *) image->imageData ) + image->roi->yOffset * width_step + image->roi->xOffset;
			}
			else
			{
				image_data = ( (const char *) image->imageData );
			}
			return get_image_mean ( image_data,
									mask_data,
									width, 
									height, 
									width_step,
									mask_width_step,
									mask_value );
		}
		break;
		case ( IPL_DEPTH_16U ):
		{
			const unsigned short int * image_data;
			width_step = image->widthStep / sizeof( short int );
			if ( image->roi != NULL )
			{
				image_data = ( (const unsigned short int *) image->imageData ) + image->roi->yOffset * width_step + image->roi->xOffset;
			}
			else
			{
				image_data = ( (const unsigned short int *) image->imageData );
			}
			return get_image_mean ( image_data,
									mask_data,
									width, 
									height, 
									width_step,
									mask_width_step,
									mask_value );
		}
		break;
		case ( IPL_DEPTH_16S ):
		{
			const short int * image_data;
			width_step = image->widthStep / sizeof( short int );
			if ( image->roi != NULL )
			{
				image_data = ( (const short int *) image->imageData ) + image->roi->yOffset * width_step + image->roi->xOffset;
			}
			else
			{
				image_data = ( (const short int *) image->imageData );
			}
			return get_image_mean ( image_data,
									mask_data,
									width, 
									height, 
									width_step,
									mask_width_step,
									mask_value );
		}
		break;
		case ( IPL_DEPTH_32S ):
		{
			const long int * image_data;
			width_step = image->widthStep / sizeof( long int );
			if ( image->roi != NULL )
			{
				image_data = ( (const long int *) image->imageData ) + image->roi->yOffset * width_step + image->roi->xOffset;
			}
			else
			{
				image_data = ( (const long int *) image->imageData );
			}
			return get_image_mean ( image_data,
									mask_data,
									width, 
									height, 
									width_step,
									mask_width_step,
									mask_value );
		}
		break;
		case ( IPL_DEPTH_32F ):
		{
			const float * image_data;
			width_step = image->widthStep / sizeof( float );
			if ( image->roi != NULL )
			{
				image_data = ( (const float *) image->imageData ) + image->roi->yOffset * width_step + image->roi->xOffset;
			}
			else
			{
				image_data = ( (const float *) image->imageData );
			}
			return get_image_mean ( image_data,
									mask_data,
									width, 
									height, 
									width_step,
									mask_width_step,
									mask_value );
		}
		break;
		case ( IPL_DEPTH_64F ):
		{
			const double * image_data;
			width_step = image->widthStep / sizeof( double );
			if ( image->roi != NULL )
			{
				image_data = ( (const double *) image->imageData ) + image->roi->yOffset * width_step + image->roi->xOffset;
			}
			else
			{
				image_data = ( (const double *) image->imageData );
			}
			return get_image_mean ( image_data,
									mask_data,
									width, 
									height, 
									width_step,
									mask_width_step,
									mask_value );
		}
		break;
		default:
			return 0;
		break;
	}
}

template<class type> double get_image_variance ( 	const type * image_data,
													const unsigned char * mask_data,
													unsigned int width,
													unsigned int height,
													unsigned int width_step,
													unsigned int mask_width_step,
													unsigned char mask_value,
													double mean )
{
	double 	sum = 0,
			tmp;
	unsigned int nb_pixels = 0; 
	for ( unsigned int i = 0; i < height; ++ i )
	{
		for ( unsigned int j = 0; j < width; ++ j )
		{
			if ( mask_data[ i * mask_width_step + j] == mask_value )
			{
				++ nb_pixels;
				tmp = image_data[ i * width_step + j ] - mean;
				sum += tmp * tmp;
			}
		}
	}
	if ( nb_pixels == 0 )
		return nan("");
	sum /= ( nb_pixels );
	return sum;
	
}

GET_IMAGE_VARIANCE_TEMP_2( unsigned char )
GET_IMAGE_VARIANCE_TEMP_2( char )

GET_IMAGE_VARIANCE_TEMP_2( unsigned short int )
GET_IMAGE_VARIANCE_TEMP_2( short int )

GET_IMAGE_VARIANCE_TEMP_2( unsigned long int )
GET_IMAGE_VARIANCE_TEMP_2( long int )

GET_IMAGE_VARIANCE_TEMP_2( float )
GET_IMAGE_VARIANCE_TEMP_2( double )

double get_image_variance( 	const IplImage * image,
							const IplImage * mask,
							unsigned char mask_value,
							double mean  )
{
unsigned int width,
				 width_bis,
				 height_bis,
				 height;
	unsigned int width_step;
	if ( mask->depth != IPL_DEPTH_8U )
		return 0;
	
	
	if ( image->roi != NULL )
	{
		width = image->roi->width;
		height = image->roi->height;
	}
	else
	{
		width = image->width;
		height = image->height;
	}
	
	const unsigned char * mask_data;
	unsigned int mask_width_step = mask->widthStep;
	if ( mask->roi != NULL )
	{
		width_bis = mask->roi->width;
		height_bis = mask->roi->height;
		mask_data = ( const unsigned char * ) mask->imageData + mask->roi->yOffset * mask_width_step + mask->roi->xOffset;
	}
	else
	{
		width_bis = mask->width;
		height_bis = mask->height;
		mask_data = ( const unsigned char * ) mask->imageData ;
	}
	
	if ( width_bis != width || height_bis != height )
		return 0;
	
	switch( image->depth )
	{
		case ( IPL_DEPTH_8U ):
		{
			const unsigned char * image_data;
			width_step = image->widthStep;
			if ( image->roi != NULL )
			{
				image_data = ( (const unsigned char *) image->imageData ) + image->roi->yOffset * width_step + image->roi->xOffset;

			}
			else
			{
				image_data = ( (const unsigned char *) image->imageData );
			}
			return get_image_variance (	image_data,
										mask_data,
										width, 
										height, 
										width_step,
										mask_width_step,
										mask_value,
										mean );
		}
		break;
		case ( IPL_DEPTH_8S ):
		{
			const char * image_data;
			width_step = image->widthStep;
			if ( image->roi != NULL )
			{
				image_data = ( (const char *) image->imageData ) + image->roi->yOffset * width_step + image->roi->xOffset;
			}
			else
			{
				image_data = ( (const char *) image->imageData );
			}
			return get_image_variance (	image_data,
										mask_data,
										width, 
										height, 
										width_step,
										mask_width_step,
										mask_value,
										mean );
		}
		break;
		case ( IPL_DEPTH_16U ):
		{
			const unsigned short int * image_data;
			width_step = image->widthStep / sizeof( short int );
			if ( image->roi != NULL )
			{
				image_data = ( (const unsigned short int *) image->imageData ) + image->roi->yOffset * width_step + image->roi->xOffset;
			}
			else
			{
				image_data = ( (const unsigned short int *) image->imageData );
			}
			return get_image_variance (	image_data,
										mask_data,
										width, 
										height, 
										width_step,
										mask_width_step,
										mask_value,
										mean );
		}
		break;
		case ( IPL_DEPTH_16S ):
		{
			const short int * image_data;
			width_step = image->widthStep / sizeof( short int );
			if ( image->roi != NULL )
			{
				image_data = ( (const short int *) image->imageData ) + image->roi->yOffset * width_step + image->roi->xOffset;
			}
			else
			{
				image_data = ( (const short int *) image->imageData );
			}
			return get_image_variance (	image_data,
										mask_data,
										width, 
										height, 
										width_step,
										mask_width_step,
										mask_value,
										mean );
		}
		break;
		case ( IPL_DEPTH_32S ):
		{
			const long int * image_data;
			width_step = image->widthStep / sizeof( long int );
			if ( image->roi != NULL )
			{
				image_data = ( (const long int *) image->imageData ) + image->roi->yOffset * width_step + image->roi->xOffset;
			}
			else
			{
				image_data = ( (const long int *) image->imageData );
			}
			return get_image_variance (	image_data,
										mask_data,
										width, 
										height, 
										width_step,
										mask_width_step,
										mask_value,
										mean );
		}
		break;
		case ( IPL_DEPTH_32F ):
		{
			const float * image_data;
			width_step = image->widthStep / sizeof( float );
			if ( image->roi != NULL )
			{
				image_data = ( (const float *) image->imageData ) + image->roi->yOffset * width_step + image->roi->xOffset;
			}
			else
			{
				image_data = ( (const float *) image->imageData );
			}
			return get_image_variance (	image_data,
										mask_data,
										width, 
										height, 
										width_step,
										mask_width_step,
										mask_value,
										mean );
		}
		break;
		case ( IPL_DEPTH_64F ):
		{
			const double * image_data;
			width_step = image->widthStep / sizeof( double );
			if ( image->roi != NULL )
			{
				image_data = ( (const double *) image->imageData ) + image->roi->yOffset * width_step + image->roi->xOffset;
			}
			else
			{
				image_data = ( (const double *) image->imageData );
			}
			return get_image_variance (	image_data,
										mask_data,
										width, 
										height, 
										width_step,
										mask_width_step,
										mask_value,
										mean );
		}
		break;
		default:
			return 0;
		break;
	}
}





