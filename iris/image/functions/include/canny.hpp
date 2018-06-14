
#ifndef _CANNY_HPP_
	#define _CANNY_HPP_
	#include <opencv/cv.h>
	/**@fn
	 * @param grad : amplitude du gradient
	 * @param image_in : image
	 * @param sigma: variance du noyau gaussien
	 * @param scaling : coefficient de réduction de l'image
	 * @brief
	 * Filtre de Canny.
	 * 
	 **/
	int canny ( IplImage * grad,
				IplImage * ori,
				const IplImage * image_in,
				double sigma,
				double scaling,
				double vert,
				double horz );

	
#endif
