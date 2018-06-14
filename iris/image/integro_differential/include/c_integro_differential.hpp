
#ifndef _C_INTEGRO_DIFFERENTIAL_HPP_
	#define _C_INTEGRO_DIFFERENTIAL_HPP_
	#include <opencv/cv.h>
	#include "../../polar/lib_polar.hpp"
	#define C_INTEGRODIFF_SEARCH_CIRCLE(type) \
	template double c_integro_differential_operator :: search_best_circle_params (	const type * data,\
																							const unsigned char * mask,\
																							unsigned int width,\
																							unsigned int height,\
																							unsigned int width_step,\
																							unsigned int mask_width_step,\
																							double x_min,\
																							double x_max,\
																							unsigned int n_x,\
																							double y_min,\
																							double y_max,\
																							unsigned int n_y,\
																							double r_min,\
																							double r_max );
	
	/**@class
	 * @brief
	 * Classe qui permet d'utiliser l'opérateur intégro différentiel.
	 * 
	 **/
	class c_integro_differential_operator
	{
		public:
			c_integro_differential_operator();
		
			/**@fn 
			 * @param width : largeur maximale d'une image
			 * @param height : hauteur maximale d'une image
			 * @param nb_directions : nombre de directions maximum
			 * @param nb_samples : nombre d'échantillon dans une direction
			 * @param kernel_size : taille du noyau gaussian
			 * 
			 **/
			c_integro_differential_operator ( unsigned int width,
											  unsigned int height,
											  unsigned int nb_directions,
											  unsigned int nb_samples,
											  unsigned int kernel_size );
			/**@fn 
			 * @param width : largeur maximale d'une image
			 * @param height : hauteur maximale d'une image
			 * @param nb_directions : nombre de directions maximum
			 * @param nb_samples : nombre d'échantillon dans une direction
			 * @param kernel_size : taille du noyau gaussian
			 * 
			 **/
			void setup (	unsigned int width,
							unsigned int height,
							unsigned int nb_directions,
							unsigned int nb_samples,
							unsigned int kernel_size );
			
			/**@fn 
			 * @param data: pixels de l'image
			 * @param mask : masque de l'image
			 * @param width : largeur de l'image
			 * @param height : hauteur de l'image
			 * @param width_step : taille réelle d'une ligne de l'image
			 * @param x_min, x_max, y_min, y_max : zone de recherche
			 * @param r_min : rayon minimal
			 * @param r_max : rayon maximal
			 * @brief
			 * Op. intégro différentiel
			 **/
			template <class type> double search_best_circle_params (	const type * data,
																	const unsigned char * mask,
																	unsigned int width,
																	unsigned int height,
																	unsigned int width_step,
																	unsigned int mask_width_step,
																	double x_min,
																	double x_max,
																	unsigned int n_x,
																	double y_min,
																	double y_max,
																	unsigned int n_y,
																	double r_min,
																	double r_max );
			/**@fn 
			 * @param image: pixels de l'image
			 * @param mask : masque de l'image
			 * @param x_min, x_max, y_min, y_max : zone de recherche
			 * @param r_min : rayon minimal
			 * @param r_max : rayon maximal
			 * @brief
			 * Op. intégro différentiel
			 **/
			double search_best_circle_params (	const IplImage * image,
												const IplImage * mask,
												double x_min,
												double x_max,
												unsigned int n_x,
												double y_min,
												double y_max,
												unsigned int n_y,
												double r_min,
												double r_max );
					
			
			
			/**@fn
			 * @brief
			 * Destructeur
			 * 
			 **/
			~c_integro_differential_operator();
		
		//Accesseurs
			inline double x() const
			{
				return _x;
			}
			inline double y() const
			{
				return _y;
			}
			inline double r() const
			{
				return _r;
			}
			
		protected:
			
			void compute_integral( );
			
			void compute_derivate( );
			
			unsigned int select_best_radius( );
			
			void free();
			void initialize();
			
			//Paramètres optimaux
			double	_x, 
					_y,
					_r;
					
			//Transformée polaire
			c_polar * polar_obj;
			unsigned int _nb_directions,
						 _nb_samples;
			
			IplImage * polar_image,
					 * polar_mask;
						 
				double * _polar_data;
				unsigned int _polar_image_width_step;
				
				unsigned char * _polar_mask;
				unsigned int _polar_mask_width_step;
			
			//Noyau gaussian + maximisation
			double * integrals,
				   * derivates;
			
			double * gaussian_kernel;
			unsigned int _kernel_size;
			IplConvKernel * structuring_element_1;
			
			
	};
	
	
	
	
	
	
	
#endif
