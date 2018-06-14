#include "image_2_polar.hpp"
#include <iostream>
using namespace std;
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
											unsigned char * p_mask,
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
											unsigned int mask_width_step ) //
{
	double x,
		   y;

	double 	r,
			theta = 0;

	double 	r_step,
			theta_step;
			
	//Calcul des pas
	r_step = ( r_max - r_min) / ( (double) nb_radii );
	theta_step = 2 * M_PI / nb_dir;
	
	r = r_min;
	for ( unsigned int i = 0; i < nb_radii; ++ i )
	{
		theta = 0;
		for ( unsigned int j = 0; j < nb_dir; ++ j )
		{
			x = r * cos(theta) + x_center;
			y = r * sin(theta) + y_center;
			
			interpol_2d ( 	p_data[ i * p_img_width_step + j ],
							p_mask[ i * p_mask_width_step  + j ],
							x,
							y,
							img_data,
							img_mask,
							width,
							height,
							width_step,
							mask_width_step ); // 
			theta += theta_step;
		}
		r += r_step;
	}
}

IMAGE_2_POLAR_DEF(unsigned char)
IMAGE_2_POLAR_DEF(char)
IMAGE_2_POLAR_DEF(unsigned short int)
IMAGE_2_POLAR_DEF(short int)
IMAGE_2_POLAR_DEF(unsigned long int)
IMAGE_2_POLAR_DEF(long int)
IMAGE_2_POLAR_DEF(float)
IMAGE_2_POLAR_DEF(double)















