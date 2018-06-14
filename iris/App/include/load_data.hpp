#ifndef _LOAD_DATA_HPP_
	#define _LOAD_DATA_HPP_
	#include "lib_image.hpp"
	#include <opencv/cv.h>
	#include <cstdio>
	struct i_code_data
	{
		i_code_data();
		unsigned int	nb_directions,
						nb_samples;
		int type;
		IplImage * image, * mask;
		
		
		~i_code_data();
	};
	
	
	int fread_iris_code(	FILE * file,
							struct i_code_data * data );
	
	int fread_iris_code(	FILE * file,
							struct i_code_data * data,
							int n );
	
	
#endif
