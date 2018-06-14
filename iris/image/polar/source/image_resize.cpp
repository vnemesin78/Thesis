#include "image_resize.hpp"
#include <iostream>
using namespace std;
template <class type> void image_resize ( 	type * data_out,
											unsigned char * mask_out,
											unsigned int width_out,
											unsigned int height_out,
											unsigned int width_step_out, 
											unsigned int mask_width_step_out,
											const type * data_in, 
											const unsigned char * mask_in,
											unsigned int width_in,
											unsigned int height_in,
											unsigned int width_step_in,
											unsigned int mask_width_step_in )
{
	double x, dx,
		   y, dy;
		   
	dx = ((double) width_in) / width_out;
	dy = ((double) height_in) / height_out;
	for ( unsigned int i = 0; i < height_out; ++ i )
	{
		y = i * dy;
		for ( unsigned int j = 0; j < width_out; ++ j )
		{
			x = j * dx;
			resample_2d(	data_out[ i * width_step_out + j ],
							mask_out[ i * mask_width_step_out + j ],
							x,
							y,
							dx,
							dy,
							data_in,
							mask_in,
							width_in,
							height_in,
							width_step_in,
							mask_width_step_in );

		}
	}
}

IMAGE_RESIZE_DEF(unsigned char)
IMAGE_RESIZE_DEF(char)

IMAGE_RESIZE_DEF(unsigned long int)
IMAGE_RESIZE_DEF(long int)

IMAGE_RESIZE_DEF(unsigned short int)
IMAGE_RESIZE_DEF(short int)

IMAGE_RESIZE_DEF(double)
IMAGE_RESIZE_DEF(float)







