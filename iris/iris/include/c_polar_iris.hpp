/**@file c_polar_iris.hpp
 * @brief 
 * Ce fichier contient la classe qui permet de calculer la transformée polaire normalisée d'un iris.
 * 
 **/
#ifndef _C_POLAR_IRIS_
	#define _C_POLAR_IRIS_
	#include <opencv/cv.h>
	#include <opencv/highgui.h>
	#include <cmath>
	#include "lib_image.hpp"
	#include "lib_api.hpp"
	#include <iostream>
	#include <sstream>
	using namespace std;
	
	#define C_POLAR_IRIS_COMPUTE(type) \
	template int c_polar_iris :: compute( 	const type * img_data,\
												const unsigned char * iris_mask,\
												unsigned int width,\
												unsigned int height,\
												unsigned int img_width_step,\
												unsigned int mask_width_step,\
												const double & x_p,\
												const double & y_p,\
												const double & a_p,\
												const double & b_p,\
												const double & theta_p,\
												const double & x_i,\
												const double & y_i,\
												const double & r_i );
	class c_polar_iris
	{
		public:
			/**@fn 
			 * @brief
			 * Constructeur.
			 **/
			c_polar_iris( void );
		
			/**@fn 
			 * @brief
			 * Constructeur.
			 **/
			c_polar_iris(	unsigned int nb_directions,
							unsigned int nb_samples,
							unsigned int nb_samples_iris,
							unsigned int width,
							unsigned int height,
							 ostream * _err_stream = NULL );
			
			int default_setup ( unsigned int width,
								unsigned int height );
			
			int setup ( 	api_parameters & params,
							unsigned int width,
							unsigned int height,
							 ostream * _err_stream = NULL,
							const char * n_space = "polar",
							const char * nb_dir_name = "nb_directions",
							const char * nb_samples_name = "nb_samples",
							const char * nb_samples_iris_name = "nb_samples_iris");
			
			
			
			/**@fn
			 * @brief
			 * Setup
			 */
			int setup (		unsigned int nb_directions,
							unsigned int nb_samples,
							unsigned int nb_samples_iris,
							unsigned int width,
							unsigned int height,
							 ostream * _err_stream = NULL );
	
			/**@fn 
			 * @brief
			 * Destructeur de poulets terminators!
			 **/
			~c_polar_iris();
		
		
			/**@fn
			 * @param img_data : pixels de l'image
			 * @param iris_mask : masque de l'image
			 * @param width : largeur de l'image
			 * @param height : hauteur de l'image
			 * @param width_step : taille réelle d'une ligne de l'image
			 * @param x_p, y_p, r_p : paramètres de la pupille
			 * @param x_i, y_i, r_i : paramètres de l'iris
			 * @brief
			 * Cette méthode calcule l'iris code de l'image.
			 **/
			template <class type> int compute ( 	const type * img_data,
													const unsigned char * iris_mask,
													unsigned int width,
													unsigned int height,
													unsigned int img_width_step,
													unsigned int mask_width_step,
													const double & x_p,
													const double & y_p,
													const double & a_p,
													const double & b_p,
													const double & theta_p,
													const double & x_i,
													const double & y_i, 
													const double & r_i );
								
			int compute (	const IplImage * image,
							const IplImage * mask,
							const double & x_p,
							const double & y_p,
							const double & r_p,
							const double & x_i,
							const double & y_i, 
							const double & r_i );								
									
			int compute (	const IplImage * image,
							const IplImage * mask,
							const double & x_p,
							const double & y_p,
							const double & a_p,
							const double & b_p,
							const double & theta_p,
							const double & x_i,
							const double & y_i, 
							const double & r_i );					
								
								
		
		
		
		
		//Accesseurs
			/**@fn
			 * @return
			 * Transformée polaire de l'iris.
			 * 
			 **/
			inline const IplImage * polar_image() const
			{
				return _polar_image;
			}
			/**@fn 
			 * @return 
			 * masque de l'image
			 **/
			inline IplImage * polar_mask()
			{
				return _polar_mask;
			}
			/**@fn 
			 * @return 
			 * nombre d'échantillons par rayon
			 **/
			inline unsigned int nb_samples() const
			{
				return _nb_samples;
			}
			
			/**@fn 
			 * @return 
			 * nombre d'échantillons par rayon
			 **/
			inline unsigned int nb_samples_iris() const
			{
				return _nb_samples_iris;
			}
			
			/**@fn 
			 * @return 
			 * nombre de directions
			 **/
			inline unsigned int nb_directions() const
			{
				return _nb_directions;
			}
			
			/**@fn
			 * @param error_str : flux d'erreur
			 * @brief
			 * Modifie le flux d'erreur.
			 * 
			 */
			inline void set_error_stream( ostream & error_str = cout )
			{
				err_stream = &error_str;
			}
			
			inline const double * pupil_radii() const
			{
				return _pupil_radii;
			}
			
			inline const double * iris_radii() const
			{
				return _iris_radii;
			}
			
		protected:
		
			unsigned int 	_nb_directions,
							_nb_samples_iris,
							_nb_samples;
							
		protected:
			/**@fn
			 * @brief
			 * Calcul le r(theta) pour la pupille.
			 * 
			 */
			void compute_pupil_radii( 	const double & x_p,
										const double & y_p,
										const double & a_p,
										const double & b_p,
										const double & theta_p );
										
			/**@fn
			 * @brief
			 * Calcul le r(theta) pour l'iris.
			 * 
			 */								
			void compute_iris_radii(	const double & x_i,
										const double & y_i,
										const double & r_i,										
										const double & x_p,
										const double & y_p );
		
			/**@fn
			 * @brief
			 * Calcul le r(theta) pour la transformée polaire
			 * 
			 */
			void compute_f_radii( );										
										
										
										
										
										
										
										
										
			double * _pupil_radii;
			double * _iris_radii;
			double * _f_radii;
			
			c_polar * polar_obj;
			IplImage * _tmp_polar_image,
					 * _tmp_polar_mask,
					 * _polar_image,
					 * _polar_mask;
			
			ostream * err_stream;
			
			
			
			
			
			/**@fn
			 * @brief
			 * Lib. mémoire
			**/
			void free();
			
			/**@fn
			 * @brief
			 * Ini. objet
			**/
			void initialize();
		
		
		
		

	};
	
	
#endif
