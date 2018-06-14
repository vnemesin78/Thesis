/**@file c_polar.hpp
 * 
 */
#ifndef _C_POLAR_HPP_
	#define _C_POLAR_HPP_
	#include "image_2_polar.hpp"
	#include "image_resize.hpp"
	#include <opencv/cv.h>
	#include <iostream>
	using namespace std;
	#define C_POLAR_DEF(type) template int  c_polar :: compute (	const double & x_center,\
																	const double & y_center,\
																	unsigned int nb_radii,\
																	unsigned int nb_dir,\
																	unsigned int r_min,\
																	unsigned int r_max,\
																	const type * img_data,\
																	const unsigned char * img_mask,\
																	unsigned int width,\
																	unsigned int height,\
																	unsigned int width_step,\
																	unsigned int width_step_mask );
																	
	
	
	/**@class 
	 * @brief
	 * Cette classe permet de réaliser une transformée polaire d'une image.
	 * 
	 */
	class c_polar
	{
		public:
			/**@fn
			 * @param nb_radii_max : nombre de rayon maximum
			 * @param nb_dir_max : nombre de direction maximum
			 * @brief
			 * Constructeur
			 * 
			 ***/
			c_polar ( unsigned int nb_radii_max,
					  unsigned int nb_dir_max );

			/**@fn
			 * @param nb_radii_max : nombre de rayon maximum
			 * @param nb_dir_max : nombre de direction maximum
			 * @brief
			 * Setup
			 * 
			 ***/
			void setup ( unsigned int nb_radii_max,
						 unsigned int nb_dir_max );
			/**@fn 
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
			 * 
			 **/
			template <class type> int compute (	const double & x_center,
												const double & y_center,
												unsigned int nb_radii,
												unsigned int nb_dir,
												unsigned int r_min,
												unsigned int r_max,
												const type * img_data,
												const unsigned char * img_mask,
												unsigned int width,
												unsigned int height,
												unsigned int width_step,
												unsigned int width_step_mask );
												
												
			/**@fn 
			 * @param x_center : abscisse du centre de la transformée polaire
			 * @param y_center : ordonnée du centre de la transformée polaire
			 * @param nb_radii : nombre de rayon de la transformée polaire
			 * @param nb_dir : nombre de directions
			 * @param r_min : rayon min
			 * @param r_max : rayon max
			 * @param image : Image
			 * @param mask : masque de l'image
			 * @brief
			 * Cette fonction calcule la transformée polaire de l'image.
			 * 
			 **/
			int compute ( const double & x_center,
							const double & y_center,
							unsigned int nb_radii,
							unsigned int nb_dir,
							unsigned int r_min,
							unsigned int r_max,
							const IplImage * image,
							const IplImage * mask );
												
												
												
												
			/**@fn 
			 * @param p_data : image de sortie
			 * @param nb_radii : nombre de rayon
			 * @param nb_dir : nombre de direction
			 * 
			 **/
			int resize ( double * p_data,
						 unsigned char * p_mask,
						 unsigned int nb_radii,
						 unsigned int nb_dir,
						 unsigned int width_step,
						 unsigned int mask_width_step );
			
			/**@fn
			 * @param p_image : image 
			 * @param p_mask : mask
			 * @brief
			 * Redimensionne l'image polaire
			 * 
			 */
			int resize( IplImage * p_image,
						IplImage * p_mask );
			
			
			
			
			
			
			/**@fn
			 * @brief
			 * Destructeur
			 * 
			 */
			~c_polar();
		
		//Accesseurs
			/**@fn
			 * @return
			 * Nb. directions allouées
			 * 
			 **/
			inline unsigned nb_directions_max() const
			{
				return _nb_dir_max;
			}
			
			/**@fn
			 * @return
			 * Nb. directions
			 * 
			 **/
			inline unsigned nb_directions() const
			{
				return _nb_dir;
			}
			/**@fn
			 * @return
			 * Nb. rayons alloués
			 * 
			 **/
			inline unsigned nb_radii_max() const
			{
				return _nb_radii_max;
			}
			
			/**@fn
			 * @return
			 * Nb. rayons
			 * 
			 **/
			inline unsigned nb_radii() const
			{
				return _nb_radii;
			}
			
			/**@fn
			 * @return Pixel de la transformée polaire
			 */
			inline const double * image_data() const
			{
				return data;
			}
		
			/**@fn
			 * @return Pixel du masque
			 */
			inline const unsigned char * image_mask() const
			{
				return mask;
			}
		
			/**@fn
			 * @brief
			 * Renvoie l'image polaire
			 * 
			 **/
			inline const IplImage * img() const
			{
				return polar_image;
			}
		
		
			/**@fn
			 * @brief
			 * Renvoie le mask
			 * 
			 **/
			inline const IplImage * p_mask() const
			{
				return polar_mask;
			}
		
			/**@fn 
			 * @return Taille d'une ligne du masque
			 **/
			inline unsigned int mask_width_step() const
			{
				return polar_mask->widthStep;
			}
		
			/**@fn 
			 * @return Taille d'une ligne du masque
			 **/
			inline unsigned int image_width_step() const
			{
				return polar_image->widthStep / sizeof(double) ;
			}
		
		
		protected:
			
			
			/**@fn
			 * @brief
			 * Tout à zéro
			 * 
			 **/
			void initialize();
			/**@fn
			 * @brief
			 * Lib. mémoire
			 * 
			 **/
			void free();
			
			
			IplImage * polar_image;
			IplImage * polar_mask;
				double * data;
				unsigned char * mask;
			
			
			
			
			unsigned int _nb_dir_max,
						 _nb_dir,
						 _nb_radii_max,
						 _nb_radii;
			
			ostream * err_stream;
	};
	
#endif
