#ifndef _CHECK_PIXEL_HPP_
	#define _CHECK_PIXEL_HPP_

	/**@fn 
	 * @param x : abscisse du pixel
	 * @param y : ordonnée du pixel
	 * @param img_mask : masque de l'image
	 * @param width : largeur de l'image
	 * @param height : hauteru de l'image
	 * @param width_step : largeur réelle d'une ligne de l'image
	 * @brief
	 * Cette fonction dit si un pixel de l'image est valide.
	 * 
	 * 
	 **/
	inline int check_pixel (	unsigned int x,
								unsigned int y,
								const unsigned char * img_mask,
								unsigned int width,
								unsigned int height,
								unsigned int width_step )
	{
		if ( x >= width || y >= height )
			return 0;
		return ( img_mask[ x + y * width_step ] );
	}


#endif
