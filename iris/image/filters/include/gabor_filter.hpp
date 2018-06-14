
#ifndef _GABOR_FILTER_HPP_
	#define _GABOR_FILTER_HPP_
	#include <fftw3.h>
	/**@fn void create_gabor_filter (	fftw_complex * gabor_filter,
										unsigned int width,
										double min_wavelenght,
										double sigma );
	 * @param gabor_filter : réponse fréquentielle du filtre de Gabor
	 * @param width : largeur
	 * @param height : hauteur
	 * @param min_wavelenght : bande passante du filtre de base
	 * @param sigma : paramètre du filtre de Gabor
	 * @brief
	 * Cette fonction crée la réponse impulsionnelle du filtre de Gabor.
	 **/
	void create_gabor_filter (	fftw_complex * gabor_filter,
								unsigned int width,
								double min_wavelenght,
								double sigma );
	
	
	
	
#endif
