
#ifndef _IMAGE_UTILITY_HPP_42_
	#define _IMAGE_UTILITY_HPP_42_
	#include <opencv/cv.h>
	/**@def CORR_VALUE( value, min, max)
	 * @param value : value
	 * @param min : min
	 * @param max : max
	 * @brief
	 * Correct sthe value. if ( value < min) value = min; else if ( value >= max ) value = max -1;
	 * 
	 **/
	#define CORR_VALUE( value, min, max)	\
	if ( ( value ) >= ( max ) ) \
		( value ) = ( max ) - 1;\
	else if ( ( value ) < ( min ) ) \
		( value ) = ( min );
		
	#define CORR_RECT( x_start, y_start, x_end, y_end, x_offset, y_offset, width, height ) \
	CORR_VALUE ( x_start, x_offset, x_offset + width ) \
	CORR_VALUE ( x_end, x_offset, x_offset + width ) \
	CORR_VALUE ( y_start, y_offset, y_offset + height ) \
	CORR_VALUE ( y_end, y_offset, y_offset + height )
	
		
												
	
	#define GET_IMAGE_DATA( img, w, h, width_step, x_offset, y_offset, image_data )	\
	if ( ( img )->roi )\
	{\
		( w ) = ( img )->roi->width;\
		( h ) = ( img )->roi->height;\
		( x_offset ) = ( img )->roi->xOffset;\
		( y_offset ) = ( img )->roi->yOffset;\
		( width_step ) = ( img )->widthStep;\
		( image_data ) = ( unsigned char * ) ( img )->imageData + ( y_offset * ( width_step ) ) + x_offset;\
	}\
	else\
	{\
		( w ) = ( img )->width;\
		( h ) = ( img )->height;\
		( x_offset ) = 0;\
		( y_offset ) = 0;\
		( width_step ) = ( img )->widthStep;\
		( image_data ) = ( unsigned char * ) ( img )->imageData + ( y_offset * ( width_step ) ) + x_offset;\
	}
	
	inline unsigned int GET_IMAGE_WIDTH_STEP( const IplImage * image )
	{
		switch( image->depth )
		{
			case(IPL_DEPTH_8U):
				return image->widthStep;
			break;
			case(IPL_DEPTH_8S):
				return image->widthStep;
			break;
			case(IPL_DEPTH_16U):
				return image->widthStep / 2;
			break;
			case(IPL_DEPTH_16S):
				return image->widthStep / 2;
			break;
			case(IPL_DEPTH_32S):
				return image->widthStep / 4;
			break;
			case(IPL_DEPTH_32F):
				return image->widthStep / 4;
			break;
			case(IPL_DEPTH_64F):
				return image->widthStep /8;
			break;
			default:
				return image->widthStep;
			break;			
		}
		
	}
	
	inline unsigned int GET_IMAGE_WIDTH ( const IplImage * image )
	{
		if ( image->roi )
			return image->roi->width;
		else
			return image->width;
	}
	
	inline unsigned int GET_IMAGE_HEIGHT ( const IplImage * image )
	{
		if ( image->roi )
			return image->roi->height;
		else
			return image->height;
	}
	
	inline unsigned int GET_IMAGE_X_OFFSET ( const IplImage * image )
	{
		if ( image->roi )
			return image->roi->xOffset;
		else
			return 0;
	}
	
	inline unsigned int GET_IMAGE_Y_OFFSET ( const IplImage * image )
	{
		if ( image->roi )
			return image->roi->yOffset;
		else
			return 0;
	}
	
	
	
	inline char * get_image_data ( const IplImage * image )
	{
		return image->imageData + image->widthStep * GET_IMAGE_Y_OFFSET ( image ) + GET_IMAGE_X_OFFSET ( image );
	}
	
	
	#define GET_IMAGE_PIXEL( image, row, col, type )\
		( (type*) ( (image)->imageData + ( (row) + GET_IMAGE_Y_OFFSET(image) ) * ( (image)->widthStep) ) )[ (col) + GET_IMAGE_X_OFFSET(image)]
		
	
	#define GET_IMAGE_DIM( img, w, h, width_step, x_offset, y_offset )	\
	if ( ( img )->roi )\
	{\
		( w ) = ( img )->roi->width;\
		( h ) = ( img )->roi->height;\
		( x_offset ) = ( img )->roi->xOffset;\
		( y_offset ) = ( img )->roi->yOffset;\
		( width_step ) = ( img )->widthStep;\
	}\
	else\
	{\
		( w ) = ( img )->width;\
		( h ) = ( img )->height;\
		( x_offset ) = 0;\
		( y_offset ) = 0;\
		( width_step ) = ( img )->widthStep;\
	}
		
		
	/**@fn
	 * @param img_out
	 * @param img_in
	 * @brief
	 * 
	 **/
	int cv_convert_image ( 	IplImage * img_out, 
							const IplImage * img_in ); 
		
																									
#endif
