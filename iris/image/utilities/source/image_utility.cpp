#include "image_utility.hpp"

#define CONVERT_IMAGE( img_in, img_out, type1, type2)\
{\
	unsigned int 	width = GET_IMAGE_WIDTH(img_out),\
					height = GET_IMAGE_HEIGHT(img_out);\
	for ( unsigned i = 0; i < height; ++ i )\
	{\
		for ( unsigned j = 0; j < width; ++ j )\
		{\
			GET_IMAGE_PIXEL(img_out, i, j, type2) = (type2) GET_IMAGE_PIXEL(img_in, i, j, type1) ;\
		}\
	}\
}

#define CONVERT_IMAGE_2( img_in, img_out, type1)\
{\
	switch((img_out)->depth)\
	{\
		case(IPL_DEPTH_8S):\
			CONVERT_IMAGE( img_in, img_out, type1, char);\
		break;\
		case(IPL_DEPTH_8U):\
			CONVERT_IMAGE( img_in, img_out, type1, unsigned char);\
		break;\
		case(IPL_DEPTH_16S):\
			CONVERT_IMAGE( img_in, img_out, type1, short int);\
		break;\
		case(IPL_DEPTH_16U):\
			CONVERT_IMAGE( img_in, img_out, type1, unsigned short int);\
		break;\
		case(IPL_DEPTH_32S):\
			CONVERT_IMAGE( img_in, img_out, type1, long int);\
		break;\
		case(IPL_DEPTH_32F):\
			CONVERT_IMAGE( img_in, img_out, type1, float);\
		break;\
		case(IPL_DEPTH_64F):\
			CONVERT_IMAGE( img_in, img_out, type1, double);\
		break;\
		default:\
			return 1;\
		break;\
	}\
}




int cv_convert_image ( 	IplImage * img_out, 
						const IplImage * img_in )
{
	if ( img_out == NULL || img_in == NULL )
		return 1;
		
	if ( 	GET_IMAGE_WIDTH(img_in) != GET_IMAGE_WIDTH(img_out) ||
			GET_IMAGE_HEIGHT(img_in) != GET_IMAGE_HEIGHT(img_out) )
		return 1;
		

	
	switch(img_in->depth)
	{
		case(IPL_DEPTH_8S):
			CONVERT_IMAGE_2( img_in, img_out, char);
		break;
		case(IPL_DEPTH_8U):
			CONVERT_IMAGE_2( img_in, img_out, unsigned char);
		break;
		case(IPL_DEPTH_16S):
			CONVERT_IMAGE_2( img_in, img_out, short int);
		break;
		case(IPL_DEPTH_16U):
			CONVERT_IMAGE_2( img_in, img_out, unsigned short int);
		break;
		case(IPL_DEPTH_32S):
			CONVERT_IMAGE_2( img_in, img_out, long int);
		break;
		case(IPL_DEPTH_32F):
			CONVERT_IMAGE_2( img_in, img_out, float);
		break;
		case(IPL_DEPTH_64F):
			CONVERT_IMAGE_2( img_in, img_out, double);		
		break;
		default:
			return 1;
		break;
	}
	return 0;
	
	
}
