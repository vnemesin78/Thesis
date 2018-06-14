/**@file hyst_threshold.hpp
 * 
 * 
 **/
#ifndef _HYST_THRESHOLD_HPP_
	#define _HYST_THRESHOLD_HPP_
	#include <opencv/cv.h>
	/**@fn
	 * @param mask : masque
	 * @param image : image
	 * @param t1, t2 : seuils pour le filtre
	 * @brief
	 * Effectue un seuillage par histéris. La 8-connexité est utilisée.
	 * 
	 * 
	 **/
	int hyst_threshold( unsigned char * mask_data,
						unsigned int width_step_mask,
						const double * image_data,
						unsigned int width_step,
						unsigned int width,
						unsigned int height,
						double t1,
						double t2 );
	
	/**@fn
	 * @param mask : masque
	 * @param image : image
	 * @param t1, t2 : seuils pour le filtre
	 * @brief
	 * Effectue un seuillage par histéris. La 8-connexité est utilisée.
	 * 
	 * 
	 **/
	int hyst_threshold( IplImage * mask,
						const IplImage * image,
						double t1,
						double t2 );
	
#endif
