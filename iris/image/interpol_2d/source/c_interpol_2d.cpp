#include "c_interpol_2d.hpp"
#include <cstring>
#include "image_utility.hpp"
	
c_interpol_2d :: c_interpol_2d (	unsigned int width,
									unsigned int height,
									const float * factors )
{
	initialize();
	setup ( 	width, 
				height, 
				factors );
}

int c_interpol_2d :: setup (	unsigned int width,
								unsigned int height,
								const float * factors  )
{
	free();
	initialize();
	if ( ! width || ! height || ! factors )
	{
		*err_stream << "Error : Invalid argument(s) in int c_interpol_2d :: setup" << endl;
		return 1;
	}
	
	memcpy(	_factors,
			factors,
			sizeof ( float ) * 4 );
	
	_width = width;
	_height = height;
	
	_image = cvCreateImage ( 	cvSize( width, height ),
								IPL_DEPTH_64F,
								1 );
	
	_previous_mask = cvCreateImage ( 	cvSize( width, height ),
										IPL_DEPTH_8U,
										1 );
										
	_curr_mask = cvCreateImage ( 	cvSize( width, height ),
									IPL_DEPTH_8U,
									1 );
									
	structuring_element_1 = cvCreateStructuringElementEx( 3,
														  3,
														  1,
														  1,
														  CV_SHAPE_CROSS,
														  NULL);
	return 0;
}

template<class type> int c_interpol_2d :: interpolate ( 	const type * image_data,
															const IplImage * mask,
															unsigned int image_width_step,
															unsigned int nb_iter )
{
	unsigned int 	x_offset,
					y_offset,
					width,
					height,
					mask_width_step;
					
	GET_IMAGE_DIM( 	mask, 
					width, 
					height, 
					mask_width_step, 
					x_offset, 
					y_offset );
	
	if ( mask->depth != IPL_DEPTH_8U )
	{
		*err_stream << "Error : invalid mask format in int c_interpol_2d :: interpolate !" << endl;
		return 1;
	}
	if ( width > _width || height > _height )
	{
		*err_stream << "Error : invalid dimensions in int c_interpol_2d :: interpolate !" << endl;
		*err_stream << "W:" << width << "\t" << _width << endl;
		*err_stream << "H:" << height << "\t" << _height << endl;
		return 2;
	}
	
	//ROI
	cvSetImageROI( 	_curr_mask,
					cvRect(0, 0, width, height ) );
	cvSetImageROI( 	_previous_mask,
					cvRect(0, 0, width, height ) );
	cvSetImageROI( 	_image,
					cvRect(0, 0, width, height ) );
					
	//Copie du mask originale
	cvCopyImage ( mask, _curr_mask );
	
	//Copie de l'image

	memset ( 	_image->imageData,
				0,
				_image->widthStep * _image->height );
	for ( unsigned int i = 0; i < height; ++ i )
	{
		for ( unsigned int j = 0; j < width; ++ j )
		{
			( (double*) (_image->imageData + i * _image->widthStep) )[j] = (double) image_data[ i * image_width_step + j ];
		}
	}
	
	for ( unsigned int i = 0; i < nb_iter ; ++ i )
	{
		if ( interpolate( width, height ) )
			break;
	}
	
	return 0;
}
template <class type> int c_interpol_2d :: get_image ( 	type * image_data,
														unsigned int width,
														unsigned int height,
														unsigned int width_step ) const
{
	for ( unsigned int i = 0; i < height; ++ i )
	{
		for ( unsigned int j = 0; j < width; ++ j )
		{
				image_data[i * width_step + j] = (type) ( (double*) (_image->imageData + i * _image->widthStep) )[j];
			
		}
	}
	return 0;
	
}




int c_interpol_2d ::get_image( IplImage * image ) const
{
	unsigned int 	x_offset,
					y_offset,
					width,
					height,
					width_step;
					
	unsigned int 	x_offset_0,
					y_offset_0,
					width_0,
					height_0,
					width_step_0;
					
					
	GET_IMAGE_DIM( 	image, 
					width, 
					height, 
					width_step, 
					x_offset, 
					y_offset );
	
	GET_IMAGE_DIM( 	_image, 
					width_0, 
					height_0, 
					width_step_0, 
					x_offset_0, 
					y_offset_0 );
	
	if ( 	width_0 != width ||
			height_0 != height )
	{
		*err_stream << "Error : Different images dimension(s) in int c_interpol_2d :: get_image" << endl;
		return 1;
	}
	
	switch ( image->depth )
	{
		case( IPL_DEPTH_8U):
			return get_image ( 	( (unsigned char*) (image->imageData + y_offset * width_step) ) + x_offset,
								width,
								height,
								width_step / sizeof( unsigned char ) );
		break;
		case( IPL_DEPTH_8S):
			return get_image ( 	( (char*) (image->imageData + y_offset * width_step) ) + x_offset,
								width,
								height,
								width_step / sizeof( char ) );
		break;
		case( IPL_DEPTH_16U):
			return get_image ( 	( (unsigned short int*) (image->imageData + y_offset * width_step) ) + x_offset,
								width,
								height,
								width_step / sizeof( unsigned short int ) );
		break;
		case( IPL_DEPTH_16S):
			return get_image ( 	( (short int*) (image->imageData + y_offset * width_step) ) + x_offset,
								width,
								height,
								width_step / sizeof( short int ) );
		break;
		case( IPL_DEPTH_32S):
			return get_image ( 	( (long int*) (image->imageData + y_offset * width_step) ) + x_offset,
								width,
								height,
								width_step / sizeof( long int ) );
		break;
		case( IPL_DEPTH_32F):
			return get_image ( 	( (float*) (image->imageData + y_offset * width_step) ) + x_offset,
								width,
								height,
								width_step / sizeof( float ) );
		break;
		break;
		case( IPL_DEPTH_64F):
			return get_image ( 	( (double*) (image->imageData + y_offset * width_step) ) + x_offset,
								width,
								height,
								width_step / sizeof( double ) );
		break;
		default:
		*err_stream << "Error : image format in int c_interpol_2d :: get_image" << endl;
		break;
	}
	return 1;
}

c_interpol_2d :: ~c_interpol_2d()
{
	free();
	initialize();
}

void c_interpol_2d :: free()
{
	if ( _image )
		cvReleaseImage ( &_image );
	if ( _curr_mask ) 
		cvReleaseImage ( &_curr_mask );
	if ( _previous_mask )
		cvReleaseImage ( &_previous_mask );
	if (structuring_element_1)
		cvReleaseStructuringElement(&structuring_element_1);
}

void c_interpol_2d :: initialize()
{
	memset ( _factors,
			 0,
			 sizeof( float ) * 4 );
			
	_image = 0;
	_previous_mask = 0;
	_curr_mask = 0;
	err_stream = &cout;
	structuring_element_1 = 0;
}

int  c_interpol_2d :: interpolate( 	unsigned int width,
									unsigned int height )
{
	//Switch
	IplImage * tmp = _curr_mask;
	_curr_mask = _previous_mask;
	_previous_mask = tmp;
	
	//Dilatation
	cvMorphologyEx(	_previous_mask,
					_curr_mask,
					NULL,
					structuring_element_1,
					CV_MOP_DILATE,
					1 );
	
	//Diff entre les deux masques
	unsigned int nb_pixels_diff = 0;
	for ( unsigned int i = 0; i < height; ++ i )
	{
		for ( unsigned int j = 0; j < width; ++ j )
		{
			int q = ( (int) ( (unsigned char *) _curr_mask->imageData )[ i * _curr_mask->widthStep + j ] ) - ( (int) ( (unsigned char *) _previous_mask->imageData )[ i * _previous_mask->widthStep + j ] );
			if ( q > 0 )
			{
				( (double*) (_image->imageData + i * _image->widthStep) )[j] = 0;
				double sum_fact = 0,
					  fact = 0;
				nb_pixels_diff ++;
				if ( j > 0 )
				{
					fact = _factors[0] * ( (unsigned char *) _previous_mask->imageData )[ i * _previous_mask->widthStep + (j - 1) ];
					sum_fact += fact;
					( (double*) (_image->imageData + i * _image->widthStep) )[j] += fact * ( (double*) (_image->imageData + i * _image->widthStep) )[j - 1];
				}
				if ( j < width - 1)
				{
					fact = _factors[3] * ( (unsigned char *) _previous_mask->imageData )[ i * _previous_mask->widthStep + (j + 1) ];
					sum_fact += fact;
					( (double*) (_image->imageData + i * _image->widthStep) )[j] += fact * ( (double*) (_image->imageData + i * _image->widthStep) )[j + 1];
				}

				if ( i > 0 )
				{
					fact = _factors[1] * ( (unsigned char *) _previous_mask->imageData )[ ( i - 1 ) * _previous_mask->widthStep + (j) ];
					sum_fact += fact;
					( (double*) (_image->imageData + i * _image->widthStep) )[j] += fact * ( (double*) (_image->imageData + ( i - 1 ) * _image->widthStep) )[j];
				}
				if ( i < height - 1)
				{
					fact = _factors[2] * ( (unsigned char *) _previous_mask->imageData )[ ( i + 1 ) * _previous_mask->widthStep + (j) ];
					sum_fact += fact;
					( (double*) (_image->imageData + i * _image->widthStep) )[j] += fact * ( (double*) (_image->imageData + ( i + 1 ) * _image->widthStep) )[j];
				}
				if ( sum_fact > 0)
					( (double*) (_image->imageData + i * _image->widthStep) )[j] /= sum_fact;
				else
				{
					( (unsigned char *) _curr_mask->imageData )[ i * _curr_mask->widthStep + j ] = 0;
				}
			}
		}
	}
	
	return (nb_pixels_diff == 0 );
}

int c_interpol_2d :: interpolate( 	const IplImage * image,
									const IplImage * mask,
									unsigned int nb_iter )
{
	unsigned int 	mask_x_offset,
					mask_y_offset,
					mask_width,
					mask_height,
					mask_width_step;
					
	unsigned int 	x_offset,
					y_offset,
					image_width,
					image_height,
					image_width_step;
					
	GET_IMAGE_DIM( 	mask, 
					mask_width, 
					mask_height, 
					mask_width_step, 
					mask_x_offset, 
					mask_y_offset );
					
	GET_IMAGE_DIM( 	image, 
					image_width, 
					image_height, 
					image_width_step, 
					x_offset, 
					y_offset );
					
	if ( 	image_width != mask_width || 
			mask_height != image_height )
	{
		*err_stream << "Error : Different image and mask dimensions in int c_interpol_2d :: interpolate!" << endl;
		return 1;
	}
	
	switch ( image->depth )
	{
		case(IPL_DEPTH_8S):
			return ( interpolate ( 	( (char * ) ( image->imageData + y_offset * image_width_step ) ) + x_offset,
									mask,
									image_width_step / sizeof(char),
									nb_iter ) );
		break;
		case(IPL_DEPTH_8U):
			return ( interpolate ( 	( (unsigned char * ) ( image->imageData + y_offset * image_width_step ) ) + x_offset,
									mask,
									image_width_step / sizeof(unsigned char),
									nb_iter ) );
		break;
		case(IPL_DEPTH_16S):
			return ( interpolate ( 	( (short int * ) ( image->imageData + y_offset * image_width_step ) ) + x_offset,
									mask,
									image_width_step / sizeof(short int),
									nb_iter ) );
		break;
		case(IPL_DEPTH_16U):
			return ( interpolate ( 	( (unsigned short int * ) ( image->imageData + y_offset * image_width_step ) ) + x_offset,
									mask,
									image_width_step / sizeof(unsigned short int),
									nb_iter ) );
		break;
		case(IPL_DEPTH_32S):
			return ( interpolate ( 	( (long int * ) ( image->imageData + y_offset * image_width_step ) ) + x_offset,
									mask,
									image_width_step / sizeof(long int),
									nb_iter ) );
		break;
		case(IPL_DEPTH_32F):
			return ( interpolate ( 	( (float * ) ( image->imageData + y_offset * image_width_step ) ) + x_offset,
									mask,
									image_width_step / sizeof(float),
									nb_iter ) );
		break;
		case(IPL_DEPTH_64F):
			return ( interpolate ( 	( (double * ) ( image->imageData + y_offset * image_width_step ) ) + x_offset,
									mask,
									image_width_step / sizeof(double),
									nb_iter ) );
		break;
		default:
			*err_stream << "Error : Unsupported image format in int c_interpol_2d :: interpolate" << endl;
			return 1;
		break;
	}
}
C_INTERPOL_2D_INTERPOLATE(unsigned char)
C_INTERPOL_2D_INTERPOLATE(char)
C_INTERPOL_2D_INTERPOLATE(unsigned short int)
C_INTERPOL_2D_INTERPOLATE(short int)
C_INTERPOL_2D_INTERPOLATE(long int)
C_INTERPOL_2D_INTERPOLATE(float)
C_INTERPOL_2D_INTERPOLATE(double)

C_INTERPOL_2D_GET_IMAGE(unsigned char)
C_INTERPOL_2D_GET_IMAGE(char)
C_INTERPOL_2D_GET_IMAGE(unsigned short int)
C_INTERPOL_2D_GET_IMAGE(short int)
C_INTERPOL_2D_GET_IMAGE(long int)
C_INTERPOL_2D_GET_IMAGE(float)
C_INTERPOL_2D_GET_IMAGE(double)








