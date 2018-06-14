/**@file interpol_2d.hpp
 * @author V. Némesin (GSM team)
 * @brief
 * Ce fichier contient le prototype de la fonction qui permet d'effectuer une approximation bilinéaire d'une image.
 * 
 */
#ifndef _INTERPOL_2D_HPP_
	#define _INTERPOL_2D_HPP_
	#include "check_pixel.hpp"
	
	
	#define INTERPOL_2D_DEF(type) template void interpol_2d ( 	double & p_value,\
																unsigned char  & p_mask,\
																const double & x,\
																const double & y,\
																const type * img_data,\
																const unsigned char * img_mask,\
																unsigned int width,\
																unsigned int height,\
																unsigned int width_step,\
																unsigned int mask_width_step );
	
	
	/**@fn 
	 * @param[out] p_value : valeur du pixel
	 * @param[out] p_mask : masque de validité du pixel
	 * @param[in] x : abscisse du pixel
	 * @param[in] y : ordonnée du pixel
	 * @param[in] img_data : données de l'image
	 * @param[in] img_mask : masque de l'image
	 * @param[in] width : largeur de l'image
	 * @param[in] height : hauteur de l'image
	 * @brief
	 * Cette fonction effectue un filtrage billinéaire.
	 * 
	 */
	template <class type> void interpol_2d ( 	double & p_value,
												unsigned char  & p_mask,
												const double & x,
												const double & y,
												const type * img_data,
												const unsigned char * img_mask,
												unsigned int width,
												unsigned int height,
												unsigned int width_step,
												unsigned int mask_width_step );
		
	
	
	
#endif
