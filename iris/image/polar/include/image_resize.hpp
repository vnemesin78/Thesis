
#ifndef _IMAGE_RESIZE_HPP_
	#define _IMAGE_RESIZE_HPP_
	#include "resample_2d.hpp"
	#define IMAGE_RESIZE_DEF(type) template void image_resize ( type * data_out,\
																unsigned char * mask_out,\
																unsigned int width_out,\
																unsigned int height_out,\
																unsigned int width_step_out,\
																unsigned int mask_width_step_out,\
																const type * data_in,\
																const unsigned char * mask_in,\
																unsigned int width_in,\
																unsigned int height_in,\
																unsigned int width_step_in,\
																unsigned int mask_width_step_in );
	
	/**@fn
	 * @param data_out : pixels de l'image redimentionnée
	 * @param mask_out : masque de l'image redimentionnée
	 * @param width_out : largeur de l'image redimentionnée
	 * @param height_out : hauteur de l'image redimentionnée
	 * @param width_step_out : taille réelle d'une ligne de l'image redimentionnée
	 * @param data_in : pixels de l'image
	 * @param mask_in : masque de l'image
	 * @param width_in : largeur de l'image
	 * @param height_in : hauteur de l'image
	 * @param width_step_in : taille réelle d'une ligne de l'image
	 * @brief
	 * Cette fonction effectue un redimentionnement de l'image.
	 * 
	 **/
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
												unsigned int mask_width_step_in );
	
#endif
