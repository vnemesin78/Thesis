/**@file image_2_polar.hpp
 * @author V. Némesin (GSM team)
 * @brief
 * Transformée polaire des images
 * 
 * 
 */
#ifndef _IMAGE_2_POLAR_HPP_
	#define _IMAGE_2_POLAR_HPP_
	#include <cmath>
	#include "interpol_2d.hpp"
	#define IMAGE_2_POLAR_DEF(type) template void image_2_polar ( 	double * p_data,\
																	unsigned char  * p_mask,\
																	const double & x_center,\
																	const double & y_center,\
																	unsigned int nb_radii,\
																	unsigned int nb_dir,\
																	unsigned int p_img_width_step,\
																	unsigned int p_mask_width_step,\
																	unsigned int r_min,\
																	unsigned int r_max,\
																	const type * img_data,\
																	const unsigned char * img_mask,\
																	unsigned int width,\
																	unsigned int height,\
																	unsigned int width_step,\
																	unsigned int mask_width_step);
	
	/**@fn
	 * @param p_data : pixels de la transformée polaire de l'image
	 * @param p_mask : masque de la transformée polaire de l'image
	 * @param x_center : abscisse du centre de la transformée polaire
	 * @param y_center : ordonnée du centre de la transformée polaire
	 * @param nb_radii : nombre de rayon de la transformée polaire
	 * @param nb_dir : nombre de directions
	 * @param r_min : rayon min
	 * @param r_max : rayon max
	 * @param img_data : pixels de l'image
	 * @param img_mask : masque de l'image
	 * @param width : largeur de l'image
	 * @param height : hauteur de l'image
	 * @param width_step : taille réelle d'une ligne de l'image.
	 * @brief
	 * Cette fonction calcule la transformée polaire de l'image.
	 **/
	template <class type> void image_2_polar ( 	double * p_data,
												unsigned char  * p_mask,
												const double & x_center,
												const double & y_center, 
												unsigned int nb_radii,
												unsigned int nb_dir,
												unsigned int p_img_width_step, //
												unsigned int p_mask_width_step, //
												unsigned int r_min,
												unsigned int r_max, //inclu
												const type * img_data, 
												const unsigned char * img_mask,
												unsigned int width,
												unsigned int height,
												unsigned int width_step,
												unsigned int mask_width_step  );
#endif
