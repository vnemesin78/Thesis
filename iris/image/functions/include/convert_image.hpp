#ifndef _CONVERT_IMAGE_HPP_
	#define _CONVERT_IMAGE_HPP_

	#include <opencv/cv.h>
	#include <opencv/highgui.h>

	/**@fn
	 * @param[out] img_out : out image 
	 * @param[in] img_in : in image
	 * @brief
	 * Conversion
	 * 
	 **/
	int convert_image ( IplImage * img_out,
						const IplImage * img_in );
					
					

	template <class type1, class type2> void convert_image_base ( 	type2 * data_out,
																	const type1 * data_in,
																	unsigned int width,
																	unsigned int height,
																	unsigned int width_step_out,
																	unsigned int width_step_int );
#endif
