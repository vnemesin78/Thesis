#include "convert_image.hpp"
#include "image_utility.hpp"
	
#define _M_convert_image(type1,type2)\
template void convert_image_base ( 	type2 * data_out,\
								const type1 * data_in,\
								unsigned int width,\
								unsigned int height,\
								unsigned int width_step_out,\
								unsigned int width_step_int );
#define _M_convert_image_2(type2)\
template int convert_image_f ( 	type2 * img_out_data,\
								const IplImage * img_int,\
								unsigned int width_step_out );

#define _M_CONVERT(type)\
	_M_convert_image_2(type)\
	_M_convert_image(type, unsigned char)\
	_M_convert_image(type, char)\
	_M_convert_image(type, short int)\
	_M_convert_image(type, unsigned short int)\
	_M_convert_image(type, float)\
	_M_convert_image(type, int)\
	_M_convert_image(type, double)


template <class type1, class type2> void convert_image_base ( type2 * data_out,
														const type1 * data_in,
														unsigned int width,
														unsigned int height,
														unsigned int width_step_out,
														unsigned int width_step_in )
{
	for ( unsigned int i = 0; i < height; ++ i )
		for ( unsigned int j = 0; j < width; ++ j )
			data_out[ i * width_step_out + j] = (type2) data_in[ i * width_step_in + j];
}

template <class type2> int convert_image_f ( 	type2 * img_out_data,
												const IplImage * img_in,
												unsigned int width_step_out )
{
	
	unsigned int 	width, 
					height,
					width_step_in,
					x_offset_in,
					y_offset_in;
					
					
	GET_IMAGE_DIM( 	img_in, 
					width, 
					height, 
					width_step_in, 
					x_offset_in, 
					y_offset_in );
					
	switch (img_in->depth )
	{
		case( IPL_DEPTH_8U ):
		{
			unsigned char * img_in_data = ( (unsigned char*) ( img_in->imageData + width_step_in * y_offset_in ) ) + x_offset_in;
			width_step_in /= sizeof(char);
			convert_image_base (img_out_data, img_in_data, width, height, width_step_out, width_step_in );
		}
		break;
		case( IPL_DEPTH_8S ):
		{
			char * img_in_data = ( (char*) ( img_in->imageData + width_step_in * y_offset_in ) ) + x_offset_in;
			width_step_in /= sizeof(char);
			convert_image_base  (img_out_data, img_in_data, width, height, width_step_out, width_step_in );
		}
		break;
		case( IPL_DEPTH_16U ):
		{
			unsigned short int * img_in_data = ( (unsigned short int *) ( img_in->imageData + width_step_in * y_offset_in ) ) + x_offset_in;
			width_step_in /= sizeof(short int);
			convert_image_base  (img_out_data, img_in_data, width, height, width_step_out, width_step_in );
		}
		break;
		case( IPL_DEPTH_16S ):
		{
			short int * img_in_data = ( (short int *) ( img_in->imageData + width_step_in * y_offset_in ) ) + x_offset_in;
			width_step_in /= sizeof(short int);
			convert_image_base  (img_out_data, img_in_data, width, height, width_step_out, width_step_in );
		}
		break;
		case( IPL_DEPTH_32F ):
		{
			float * img_in_data = ( (float *) ( img_in->imageData + width_step_in * y_offset_in ) ) + x_offset_in;
			width_step_in /= sizeof(float);
			convert_image_base  (img_out_data, img_in_data, width, height, width_step_out, width_step_in );
		}
		break;
		case( IPL_DEPTH_32S ):
		{
			int * img_in_data = ( (int *) ( img_in->imageData + width_step_in * y_offset_in ) ) + x_offset_in;
			width_step_in /= sizeof(int);
			convert_image_base (img_out_data, img_in_data, width, height, width_step_out, width_step_in );
		}
		break;
		case( IPL_DEPTH_64F ):
		{
			double * img_in_data = ( (double *) ( img_in->imageData + width_step_in * y_offset_in ) ) + x_offset_in;
			width_step_in /= sizeof(double);
			convert_image_base (img_out_data, img_in_data, width, height, width_step_out, width_step_in );
		}
		break;
		default:
			return -1;
		break;
	}
	return 0;
}



int convert_image ( IplImage * img_out,
					const IplImage * img_in )
{
	unsigned int 	width_out, 
					height_out,
					width_step_out,
					x_offset_out,
					y_offset_out;
					
	unsigned int 	width_in, 
					height_in;
					
	GET_IMAGE_DIM( 	img_out, 
					width_out, 
					height_out, 
					width_step_out, 
					x_offset_out, 
					y_offset_out );
					
	width_in = GET_IMAGE_WIDTH(img_in);
	height_in = GET_IMAGE_HEIGHT(img_in);
					
	if ( width_in != width_out || height_in != height_out )
		return -1;
	
	switch (img_out->depth )
	{
		case( IPL_DEPTH_8U ):
		{
			unsigned char * img_out_data = ( (unsigned char*) ( img_out->imageData + width_step_out * y_offset_out ) ) + x_offset_out;
			width_step_out /= sizeof(char);
			return convert_image_f (img_out_data, img_in, width_step_out );
		}
		break;
		case( IPL_DEPTH_8S ):
		{
			char * img_out_data = ( (char*) ( img_out->imageData + width_step_out * y_offset_out ) ) + x_offset_out;
			width_step_out /= sizeof(char);
			return convert_image_f (img_out_data, img_in, width_step_out );
		}
		break;
		case( IPL_DEPTH_16U ):
		{
			unsigned short int * img_out_data = ( (unsigned short int *) ( img_out->imageData + width_step_out * y_offset_out ) ) + x_offset_out;
			width_step_out /= sizeof(short int);
			return convert_image_f (img_out_data, img_in, width_step_out );
		}
		break;
		case( IPL_DEPTH_16S ):
		{
			short int * img_out_data = ( (short int *) ( img_out->imageData + width_step_out * y_offset_out ) ) + x_offset_out;
			width_step_out /= sizeof(short int);
			return convert_image_f (img_out_data, img_in, width_step_out );
		}
		break;
		case( IPL_DEPTH_32F ):
		{
			float * img_out_data = ( (float *) ( img_out->imageData + width_step_out * y_offset_out ) ) + x_offset_out;
			width_step_out /= sizeof(float);
			return convert_image_f (img_out_data, img_in, width_step_out );
		}
		break;
		case( IPL_DEPTH_32S ):
		{
			int * img_out_data = ( (int *) ( img_out->imageData + width_step_out * y_offset_out ) ) + x_offset_out;
			width_step_out /= sizeof(int);
			return convert_image_f (img_out_data, img_in, width_step_out );
		}
		break;
		case( IPL_DEPTH_64F ):
		{
			double * img_out_data = ( (double *) ( img_out->imageData + width_step_out * y_offset_out ) ) + x_offset_out;
			width_step_out /= sizeof(double);
			return convert_image_f (img_out_data, img_in, width_step_out );
		}
		break;
		default:
			return -1;
		break;
	}
	return 0;
}

_M_CONVERT(unsigned char)
_M_CONVERT(char)
_M_CONVERT(unsigned short int)
_M_CONVERT(short int)
_M_CONVERT(int)
_M_CONVERT(float)
_M_CONVERT(double)
