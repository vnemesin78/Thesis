/**@file c_iris_code.hpp
 * @brief
 * Ce fichier contient le prototype de la classe qui permet de calculer un iris code.
 * 
 **/
#ifndef _C_IRIS_CODE_HPP_
	#define _C_IRIS_CODE_HPP_
	#include "lib_image.hpp"
	#include "lib_api.hpp"
	#include <opencv/cv.h>
	#include <opencv/highgui.h>
	#include <fstream>
	#define C_IRIS_CODE_COMPUTE_IRIS_CODE(type)\
	template void c_iris_code :: compute_iris_code ( 	const type * img_data,\
														const unsigned char * mask,\
														unsigned int img_width_step,\
														unsigned int mask_width_step );
	
	using namespace std;
	/**@class iris_code
	 * @brief
	 * Cette classe permet de calculer un iris code à partir d'une image
	 * 
	 **/
	class c_iris_code
	{
		public:
			c_iris_code()
			{
				initialize();
			}
		
		
			/**@fn 
			 * @param width : largeur de l'image
			 * @param height : hauteur de l'image
			 * @param nb_directions : nombre
			 * @param nb_samples : nombre d'échantillons sur chq rayons
			 * @param wavelenght 
			 * @param sigma
			 * @brief
			 * Constructeur
			 **/
			c_iris_code (	unsigned int nb_directions,
							unsigned int nb_samples,
							double wavelenght,
							double sigma,
							ostream * _err_stream = NULL );
			/**@fn 
			 * @param width : largeur de l'image
			 * @param height : hauteur de l'image
			 * @param nb_directions : nombre
			 * @param nb_samples : nombre d'échantillons sur chq rayons
			 * @param wavelenght 
			 * @param sigma
			 * @brief
			 * Setup
			 **/
			int setup (	unsigned int nb_directions = 240,
						unsigned int nb_samples = 20,
						double wavelenght = 18,
						double sigma = 0.5,
						ostream * _err_stream = NULL );
		
			int default_setup( );
		
			/**@fn
			 * @param params : paramètres
			 * @brief
			 * Setup
			 * 
			 **/
			int setup(	api_parameters & params,
						ostream * _err_stream = NULL,
						const char * prefix = "iris_code",
						const char * nb_directions_name = "nb_directions",
						const char * nb_samples_name = "nb_samples",
						const char * wavelenght_name = "wavelenght",
						const char * sigma_name = "sigma" );
		
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
			int compute_iris_code (	const IplImage * image,
									const IplImage * mask );
									
									
			template <class type> void compute_iris_code ( 	const type * img_data,
															const unsigned char * mask,
															unsigned int img_width_step,
															unsigned int mask_width_step );
		
			/**@fn 
			 * @brief
			 * Destructeur
			 * 
			 */
			~c_iris_code();
			
		//Accesseurs
			/**@fn 
			 * @return 
			 * iris code
			 **/
			inline const IplImage * code() const
			{
				return _code_img;
			}
			/**@fn 
			 * @return 
			 * masque de l'image
			 **/
			inline const IplImage * mask() const
			{
				return _mask_img;
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
			
			inline double nrj_ratio() const
			{
				return _nrj_ratio;
			}
			
		protected:
			c_interpol_2d interpol_2d_obj;
		
		
		
			fftw_plan plan_in,
					  plan_out;
		
			//Data
			IplImage * _code_img,
					 * _mask_img;
			
			IplImage * _g_img_re,
					 * _g_img_im;
			
			double _nrj_ratio;
			//Parameters
			unsigned int 	_nb_directions,
							_nb_samples,
							_nb_pixels;
			double _wavelenght,
				   _sigma;
				   
			//Tmp
			fftw_complex * _array;
			fftw_complex * _f_array;
			
			fftw_complex * tmp;
			
			//Log-Gabor filter
			fftw_complex * gabor_filter;

			IplConvKernel * structuring_element_1;
			IplImage * tmp1, * tmp2;
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
			
			//errors stream
			ostream * err_stream;
	};
	
	
#endif
