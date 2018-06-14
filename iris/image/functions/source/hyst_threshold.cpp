#include "hyst_threshold.hpp"

int hyst_threshold( unsigned char * mask_data,
					unsigned int width_step_mask,
					const double * image_data,
					unsigned int width_step,
					unsigned int width,
					unsigned int height,
					double t1,
					double t2 )
{
	
	unsigned int * x_pixels,
				 * y_pixels,
				 * x_stack,
				 * y_stack;
		
	unsigned int n_pixels = 0,
				 nb_pixels;
	
	int x_tmp[8],
		y_tmp[8];
	
	
	if ( t1 < t2 || t2 < 0 )
		return -1;
	
	nb_pixels = width * height;
	x_pixels = new unsigned int[ nb_pixels ];
	y_pixels = new unsigned int[ nb_pixels ];
	
	for ( unsigned int i = 0; i < height; ++ i )
	{
		for ( unsigned int j = 0; j < width; ++ j )
		{	
			if ( image_data [ i * width_step + j ] > t1 )
			{
				x_pixels [ n_pixels ] =  j;
				y_pixels [ n_pixels ] =  i;
				n_pixels ++;
			}
		}
	}
	
	//Stack
	x_stack = new unsigned int[nb_pixels]; //A voir pour l'optimisation
	y_stack = new unsigned int[nb_pixels];
	
	memcpy( x_stack, x_pixels, sizeof(int) * n_pixels ); 
	memcpy( y_stack, y_pixels, sizeof(int) * n_pixels );
	
	
	
	//Reset the mask
	for ( unsigned int i = 0; i < height; ++ i )
	{
		memset ( mask_data + i * width_step_mask,
				 0,
				 sizeof( char ) * width );
	}
	
	//Mark points as edge
	for ( unsigned int i = 0; i < n_pixels; ++ i )
	{
		mask_data[ y_pixels [ i ] * width_step_mask + x_pixels[ i ] ] = 255; // A voir
	}
	
	//Directions
	x_tmp[0] = -1;
	y_tmp[0] = 0;
	
	x_tmp[1] = 1;
	y_tmp[1] = 0;
	
	x_tmp[2] = -1;
	y_tmp[2] = -1;
	
	x_tmp[3] = 0;
	y_tmp[3] = -1;
	
	x_tmp[4] = 1;
	y_tmp[4] = -1;
	
	x_tmp[5] = -1;
	y_tmp[5] = 1;
	
	x_tmp[6] = 0;
	y_tmp[6] = 1;
	
	
	x_tmp[7] = 1;
	y_tmp[7] = 1;
	
	int stp = n_pixels - 1;
	while ( stp >= 0 )//% While the stack is not empty
	{
		unsigned int x, y;
		x = x_stack[stp];         //% Pop next index off the stack
		y = y_stack[stp];
		-- stp;
		
		for (unsigned i = 0; i < 8; ++ i)
		{
			int x_l,
				y_l;
				
			x_l = x + x_tmp[i];
			y_l = y + y_tmp[i];
			if ( ! ( x_l < 0 || (unsigned int) x_l >= width || y_l < 0 || (unsigned int) y_l >= height ) )
			{
				unsigned int n_pixel = x_l + width_step * y_l;
				unsigned int n_pixel_mask = x_l + width_step_mask * y_l;
				if ( image_data[ n_pixel ] > t2 && ( ! mask_data[ n_pixel_mask ] ) )
				{
					++ stp;
					x_stack[stp] = x_l;
					y_stack[stp] = y_l;
					mask_data[ n_pixel_mask ] = 255;
				}
			}
		}
	}

	//Lib. mémoire
	delete[] x_pixels;
	delete[] y_pixels;
	delete[] x_stack;
	delete[] y_stack;
	
	return 0;
}
					
int hyst_threshold( IplImage * mask,
					const IplImage * image,
					double t1,
					double t2 )
{
	unsigned int mask_width,
				 mask_height,
				 mask_width_step;
				 
	unsigned int image_width,
				 image_height,
				 image_width_step;
				 
	unsigned char * mask_data;
	const double * image_data;
	
	{
		
		mask_width_step = mask->widthStep;
		if ( mask->roi )
		{
			mask_width = mask->roi->width;
			mask_height = mask->roi->height;
			mask_data = (unsigned char*) mask->imageData + mask->roi->xOffset + mask->roi->yOffset * mask_width_step;
		}
		else
		{
			mask_width = mask->width;
			mask_height = mask->height;
			mask_data = (unsigned char*) mask->imageData;
		}
		
	}
	
	{
		
		image_width_step = image->widthStep / sizeof( double );
		if ( image->roi )
		{
			image_width = image->roi->width;
			image_height = image->roi->height;
			image_data = ( (const double *) image->imageData ) + image->roi->xOffset + image->roi->yOffset * image_width_step;
		}
		else
		{
			image_width = image->width;
			image_height = image->height;
			image_data = (const double *) image->imageData;
		}
		
	}
	if ( image_width != mask_width || image_height != mask_height  || mask->depth != IPL_DEPTH_8U || image->depth != IPL_DEPTH_64F )
		return 1;
		
	return hyst_threshold(	mask_data,
							mask_width_step,
							image_data,
							image_width_step,
							image_width,
							image_height,
							t1,
							t2 );
}
