
#ifndef _RESAMPLE_2D_HPP_
	#define _RESAMPLE_2D_HPP_
	#include "check_pixel.hpp"
	#define RESAMPLE_2D_DEF(type) template void resample_2d(	type & p_value,\
																unsigned char  & p_mask,\
																const double & x,\
																const double & y,\
																const double & dx,\
																const double & dy,\
																const type * img_data,\
																const unsigned char * img_mask,\
																unsigned int width,\
																unsigned int height,\
																unsigned int width_step,\
																unsigned int mask_width_step );
	
	/**@fn
	 * @param p_value : valeur du pixel
	 * @param p_mask : masque
	 * @param x : abscisse
	 * @param y : ordonnée
	 * @param dx : largeur du pixel
	 * @param dy : hauteur du pixel
	 * @param img_data : données de l'image
	 * @param img_mask : masque de l'image
	 * @param width : largeur de l'image
	 * @param height : hauteur de l'image
	 * @param width_step : taille réelle d'une ligne de l'image
	 * @brief
	 * Ré-échantillonnage avec filtre anti-aliasing.
	 */
	template <class type> void resample_2d(	type & p_value,
											unsigned char & p_mask,
											const double & x,
											const double & y,
											const double & dx,
											const double & dy,
											const type * img_data,
											const unsigned char * img_mask,
											unsigned int width,
											unsigned int height,
											unsigned int width_step,
											unsigned int mask_width_step );
	
#endif
