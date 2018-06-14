#include "c_fusion.hpp"

c_fusion :: c_fusion ( void )
{
	initialize();
}

c_fusion :: c_fusion ( unsigned int width,
					   unsigned int height,
					   ostream * _err_stream )
{
	initialize();
	if ( setup ( width, height, _err_stream ) )
		throw ( invalid_argument("Error in c_fusion( unsigned int width, unsigned int height )!") );
}

void c_fusion :: reset( void )
{
	if ( _fuzzy_image )
		memset (	_fuzzy_image->imageData, 
					0,
					_fuzzy_image->widthStep * _fuzzy_image->height );
	
	if ( _coef_image )
		memset (	_coef_image->imageData, 
					0,
					_coef_image->widthStep * _coef_image->height );
	_sum_coef = 0;
}

int c_fusion :: setup (	unsigned int width,
						unsigned int height,
						ostream * _err_stream )
{
	free();
	initialize();
	err_stream = _err_stream;
	if ( !width || !height )
	{
		if ( err_stream )
			*err_stream << "Error : in int c_fusion::setup( unsigned int, unsigned int), width or/and height is/are 0." << endl;
		return 1;
	}
	_width = width;
	_height = height;
	
	_fuzzy_image = cvCreateImage ( 	cvSize ( width, height ),
									IPL_DEPTH_64F,
									1 );
	
	_coef_image = cvCreateImage ( 	cvSize ( width, height ),
									IPL_DEPTH_64F,
									1 );
	reset();
	
	return 0;
}

int c_fusion :: add_image (	const IplImage * image,
							const IplImage * mask,
							double score,
							int d_x )
{
	unsigned int	img_width,
					img_height,
					img_width_step,
					img_x_offset,
					img_y_offset;
	
	unsigned int	mask_width,
					mask_height,
					mask_width_step,
					mask_x_offset,
					mask_y_offset;
	
	if ( !image || !mask )
	{
		if ( err_stream )
			*err_stream << "Error : image or/and mask is/are NULL!" << endl;
		return -1;
	}
	
	
	GET_IMAGE_DIM( 	image, 
					img_width,
					img_height,
					img_width_step, 
					img_x_offset, 
					img_y_offset );
					
	GET_IMAGE_DIM( 	mask, 
					mask_width,
					mask_height,
					mask_width_step, 
					mask_x_offset, 
					mask_y_offset );
					
	if ( 	img_width 	!= 	_width 	||
			mask_width	!=	_width	||
			img_height	!=	_height	||
			mask_height	!= 	_height  )
	{
		if ( err_stream )
		{
			*err_stream << "Error : image or/and mask dimensions!" << endl;
			*err_stream << img_width << "\t" << mask_width << "\t" << _width << endl;
			*err_stream << img_height << "\t" << mask_height << "\t" << _height << endl;
		}
		return -2;
	}
	
	if ( mask->depth != IPL_DEPTH_8U )
	{
		//*err_stream << "Error : invalid mask format!" << endl;
		return -3;
	}


	switch( image->depth )
	{
		case( IPL_DEPTH_8S):
		{
			add_image (	( (char *) ( image->imageData + img_y_offset * img_width_step ) ) + img_x_offset,
						img_width_step,
						( (unsigned char *) ( mask->imageData + mask_y_offset * mask_width_step ) ) + mask_x_offset,
						mask_width_step,
						score,
						d_x );
		}
		break;
		case( IPL_DEPTH_8U):
		{
			add_image (	( (unsigned char *) ( image->imageData + img_y_offset * img_width_step ) ) + img_x_offset,
						img_width_step,
						( (unsigned char *) ( mask->imageData + mask_y_offset * mask_width_step ) ) + mask_x_offset,
						mask_width_step,
						score,
						d_x );
		}
		break;
		case( IPL_DEPTH_16S):
		{
			add_image (	( (short int *) ( image->imageData + img_y_offset * img_width_step ) ) + img_x_offset,
						img_width_step / sizeof(short int),
						( (unsigned char *) ( mask->imageData + mask_y_offset * mask_width_step ) ) + mask_x_offset,
						mask_width_step,
						score,
						d_x );
		}
		break;
		case( IPL_DEPTH_16U):
		{
			add_image (	( (unsigned short int *) ( image->imageData + img_y_offset * img_width_step ) ) + img_x_offset,
						img_width_step / sizeof(short int),
						( (unsigned char *) ( mask->imageData + mask_y_offset * mask_width_step ) ) + mask_x_offset,
						mask_width_step,
						score,
						d_x );
		}
		break;
		case( IPL_DEPTH_32S):
		{
			add_image (	( (long int *) ( image->imageData + img_y_offset * img_width_step ) ) + img_x_offset,
						img_width_step / sizeof(long int),
						( (unsigned char *) ( mask->imageData + mask_y_offset * mask_width_step ) ) + mask_x_offset,
						mask_width_step,
						score,
						d_x );
		}
		break;
		case( IPL_DEPTH_32F):
		{
			add_image (	( (float *) ( image->imageData + img_y_offset * img_width_step ) ) + img_x_offset,
						img_width_step / sizeof(float),
						( (unsigned char *) ( mask->imageData + mask_y_offset * mask_width_step ) ) + mask_x_offset,
						mask_width_step,
						score,
						d_x );
		}
		break;
		case( IPL_DEPTH_64F):
		{
			add_image (	( (double*) ( image->imageData + img_y_offset * img_width_step ) ) + img_x_offset,
						img_width_step / sizeof(double),
						( (unsigned char *) ( mask->imageData + mask_y_offset * mask_width_step ) ) + mask_x_offset,
						mask_width_step,
						score,
						d_x );
		}
		break;
		default:
			if ( err_stream )
				*err_stream << "Error : unsupported image format!" << endl;
			return -4;
		break;
	}
	return 0;
}

template <class type> void c_fusion :: add_image (	const type * img_data,
													unsigned int img_width_step,
													const unsigned char * mask_data,
													unsigned int mask_width_step,
													double score,
													int d_x )
{
	for ( unsigned int i = 0; i < _height; ++ i )
	{
		for ( unsigned int j = 0; j < _width; ++ j )
		{
			if ( (mask_data + mask_width_step * i )[(j + d_x) % _width ] ) 
			{
				( (double*) (_fuzzy_image->imageData + i * _fuzzy_image->widthStep ) )[j] += score * (img_data + img_width_step * i )[(j + d_x) % _width ];
				( (double*) (_coef_image->imageData + i * _coef_image->widthStep ) )[j] += score;
			}
			
			
		}
	}
	_sum_coef += score;
}

int c_fusion :: compute (	IplImage * image,
							IplImage * mask,
							bool binary_mask,
							double threshold  )
{
	unsigned int	img_width,
					img_height,
					img_width_step,
					img_x_offset,
					img_y_offset;
	
	unsigned int	mask_width,
					mask_height,
					mask_width_step,
					mask_x_offset,
					mask_y_offset;
	
	if ( !image || !mask )
	{
		if ( err_stream )
			*err_stream << "Error : image or/and mask is/are NULL!" << endl;
		return -1;
	}
	
	
	GET_IMAGE_DIM( 	image, 
					img_width,
					img_height,
					img_width_step, 
					img_x_offset, 
					img_y_offset );
					
	GET_IMAGE_DIM( 	mask, 
					mask_width,
					mask_height,
					mask_width_step, 
					mask_x_offset, 
					mask_y_offset );
					
	if ( 	img_width 	!= 	_width 	||
			mask_width	!=	_width	||
			img_height	!=	_height	||
			mask_height	!= 	_height  )
	{
		if ( err_stream )
			*err_stream << "Error : image or/and mask dimensions!" << endl;
		return -2;
	}
	
	if ( mask->depth != IPL_DEPTH_8U )
	{
		if ( err_stream )
			*err_stream << "Error : invalid mask format!" << endl;
		return -3;
	}
	
	switch( image->depth )
	{
		case( IPL_DEPTH_8S):
		{
			compute (	( (char *) ( image->imageData + img_y_offset * img_width_step ) ) + img_x_offset,
						img_width_step,
						( (unsigned char *) ( mask->imageData + mask_y_offset * mask_width_step ) ) + mask_x_offset,
						mask_width_step,
						binary_mask,
						threshold );
		}
		break;
		case( IPL_DEPTH_8U):
		{
			compute (	( (unsigned char *) ( image->imageData + img_y_offset * img_width_step ) ) + img_x_offset,
						img_width_step,
						( (unsigned char *) ( mask->imageData + mask_y_offset * mask_width_step ) ) + mask_x_offset,
						mask_width_step,
						binary_mask,
						threshold );
		}
		break;
		case( IPL_DEPTH_16S):
		{
			compute (	( (short int *) ( image->imageData + img_y_offset * img_width_step ) ) + img_x_offset,
						img_width_step / sizeof(short int),
						( (unsigned char *) ( mask->imageData + mask_y_offset * mask_width_step ) ) + mask_x_offset,
						mask_width_step,
						binary_mask,
						threshold );
		}
		break;
		case( IPL_DEPTH_16U):
		{
			compute (	( (unsigned short int *) ( image->imageData + img_y_offset * img_width_step ) ) + img_x_offset,
						img_width_step / sizeof(short int),
						( (unsigned char *) ( mask->imageData + mask_y_offset * mask_width_step ) ) + mask_x_offset,
						mask_width_step,
						binary_mask,
						threshold );
		}
		break;
		case( IPL_DEPTH_32S):
		{
			compute (	( (long int *) ( image->imageData + img_y_offset * img_width_step ) ) + img_x_offset,
						img_width_step / sizeof(long int),
						( (unsigned char *) ( mask->imageData + mask_y_offset * mask_width_step ) ) + mask_x_offset,
						mask_width_step,
						binary_mask,
						threshold );
		}
		break;
		case( IPL_DEPTH_32F):
		{
			compute (	( (float *) ( image->imageData + img_y_offset * img_width_step ) ) + img_x_offset,
						img_width_step / sizeof(float),
						( (unsigned char *) ( mask->imageData + mask_y_offset * mask_width_step ) ) + mask_x_offset,
						mask_width_step,
						binary_mask,
						threshold );
		}
		break;
		case( IPL_DEPTH_64F):
		{
			compute (	( (double*) ( image->imageData + img_y_offset * img_width_step ) ) + img_x_offset,
						img_width_step / sizeof(double),
						( (unsigned char *) ( mask->imageData + mask_y_offset * mask_width_step ) ) + mask_x_offset,
						mask_width_step,
						binary_mask,
						threshold );
		}
		break;
		default:
			if ( err_stream )
				*err_stream << "Error : unsupported image format!" << endl;
			return -4;
		break;
	}
	return 0;
}

template <class type> void c_fusion :: compute (	type * img_data,
													unsigned int img_width_step,
													unsigned char * mask_data,
													int mask_width_step,
													bool binary_mask,
													double threshold )
{
	for ( unsigned int i = 0; i < _height; ++ i )
	{
		for ( unsigned int j = 0; j < _width; ++ j )
		{
			if ( ( (double*) (_coef_image->imageData + i * _coef_image->widthStep ) )[j] > 0 )
			{
				(img_data + img_width_step * i )[j] =  (type) ( ( (double*) (_fuzzy_image->imageData + i * _fuzzy_image->widthStep ) )[j]/ ( (double*) (_coef_image->imageData + i * _coef_image->widthStep ) )[j] + 0.5 );

			
			
				(mask_data + mask_width_step * i )[j] = (unsigned char) ( 255.0 * ( ( ( (double*) (_coef_image->imageData + i * _coef_image->widthStep ) )[j] + _sum_coef / 2 ) / _sum_coef ) );
			}
		}
	}
	if ( binary_mask )
	{
		double t = threshold * _sum_coef;
		for ( unsigned int i = 0; i < _height; ++ i )
		{
			for ( unsigned int j = 0; j < _width; ++ j )
			{
				if ( ( ( (double*) (_coef_image->imageData + i * _coef_image->widthStep ) ) )[j] <= t )
					(mask_data + mask_width_step * i )[j] = 0;
				else
					(mask_data + mask_width_step * i )[j] = 255;
			}
		}
	}
	
}

c_fusion :: ~c_fusion( void )
{
	free();
	initialize();
}

void c_fusion :: free()
{
	if ( _fuzzy_image )
		cvReleaseImage ( &_fuzzy_image );
	if ( _coef_image )
		cvReleaseImage ( &_coef_image );	
}

void c_fusion :: initialize()
{
	_fuzzy_image = 0;
	_coef_image = 0;
	_sum_coef = 0;
	_width = 0;
	_height = 0;
	err_stream =NULL;
}


C_FUSION_TEMPLATE_ADD_IMAGE(unsigned char)
C_FUSION_TEMPLATE_ADD_IMAGE(char)
C_FUSION_TEMPLATE_ADD_IMAGE(unsigned short int)
C_FUSION_TEMPLATE_ADD_IMAGE(short int)
C_FUSION_TEMPLATE_ADD_IMAGE(long int)
C_FUSION_TEMPLATE_ADD_IMAGE(float)
C_FUSION_TEMPLATE_ADD_IMAGE(double)


C_FUSION_TEMPLATE_COMPUTE(unsigned char)
C_FUSION_TEMPLATE_COMPUTE(char)
C_FUSION_TEMPLATE_COMPUTE(unsigned short int)
C_FUSION_TEMPLATE_COMPUTE(short int)
C_FUSION_TEMPLATE_COMPUTE(long int)
C_FUSION_TEMPLATE_COMPUTE(float)
C_FUSION_TEMPLATE_COMPUTE(double)



