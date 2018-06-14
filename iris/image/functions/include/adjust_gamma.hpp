/**@file adjust_gamma.hpp
 * 
 **/
#ifndef _ADJUST_GAMMA_HPP_
	#define _ADJUST_GAMMA_HPP_

	
	/**@fn int adjust_gamma (	unsigned char * img_data,
								unsigned int width,
								unsigned int height,
								unsigned int width_step,
								double g );
	 * @param img_data : pixels de l'image
	 * @param width : largeur de l'image
	 * @param height : hauteur de l'image
	 * @param width_step : taille réelle d'une ligne de l'image
	 * @param g : gain de correciton du gamme
	 * @return
	 * - 0 si tout se déroule bien.
	 * @brief
	 * Corrige le gamma d'une image
	 */
	int adjust_gamma ( double * img_data,
					   unsigned int width,
					   unsigned int height,
					   unsigned int width_step,
					   double g );
	
	/**@fn int adjust_gamma_2 (	unsigned char * img_data,
								unsigned int width,
								unsigned int height,
								unsigned int width_step,
								double g );
	 * @param img_data : pixels de l'image
	 * @param width : largeur de l'image
	 * @param height : hauteur de l'image
	 * @param width_step : taille réelle d'une ligne de l'image
	 * @param g : gain de correciton du gamme
	 * @return
	 * - 0 si tout se déroule bien.
	 * @brief
	 * Corrige le gamma d'une image (histogramme corrigé)
	 */
	int adjust_gamma_2 ( 	double * img_data,
							unsigned int width,
							unsigned int height,
							unsigned int width_step,
							double g );
	
	
	
#endif
