#ifndef _NON_MAX_SUP_HPP_
	#define _NON_MAX_SUP_HPP_
	#include <cmath>
	#include <opencv/cv.h>
	/**@fn
	 * @param image_data : pixels de l'image
	 * @param image_width_step : taille d'une ligne de l'image
	 * @param grad_data : pixel du gradient
	 * @param grad_width_step : taille d'une ligne du gradient
	 * @param orient_data : pixels des différents angles du gradient
	 * @param orient_data_width_step : largeur d'une ligne de l'orientation
	 * @param width : largeur
	 * @param height : hauteur
	 * @param radius : rayon.
	 * @brief
	 * 
	 * 
	 * 
	 **/
	int non_max_sup ( 	double * image_data,
						unsigned int image_width_step,
						const double * grad_data,
						unsigned int grad_width_step,
						const double * orient_data,
						unsigned int orient_width_step,
						unsigned int width,
						unsigned int height,
						double radius );
	/**@fn
	 * @param image_out : image traitée
	 * @param image_grad : amplitude du gradiant
	 * @param image_orient : orientation du gradient
	 * @param radius : rayon
	 * @brief
	 * Fonction pour supprimer les non - maximum
	 **/
	int non_max_sup ( 	IplImage * image_out,
						const IplImage * image_grad,
						const IplImage * image_orient,
						double radius );
#endif
