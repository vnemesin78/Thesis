/**@file c_eyelids_segmentation.hpp
 * @author Valérian Némesin
**/
#ifndef _C_EYELIDS_SEGMENTATION_HPP_
#define _C_EYELIDS_SEGMENTATION_HPP_
#include "lib_image.hpp"
#include "lib_api.hpp"
#include <opencv/cv.h>
#include <opencv/highgui.h>
#include <sstream>
#include <iostream>

#include <exception> //Exceptions (pour ne pas faire planter le programme en cas de problèmes de mémoire ou autre)
#include <stdexcept>
using namespace std;

/**@class
 * @brief
 * Segmentation des paupières
 * 
 */	
class c_eyelids_segmentation
{
	public:
		/**@fn
		 * @brief
		 * Constructeur
		 * 
		 */
		c_eyelids_segmentation ( void );
		
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
		
		/**@fn
		 * @param nb_samples : nombre d'échantillons radials
		 * @param nb_directions : nombre de directions
		 * @param width : largeur de l'image (width * r_ratio)
		 * @param height : hauteur de l'image (height * r_ratio )
		 * @param mean_bg : moyenne du bg
		 * @param sigma_bg : écart type du bg
		 * @param nb_iter_gem : nombre d'itération de l'algo de Mouloud
		 * @param nb_iter_icm : nombre d'itération de l'algo ICM
		 * @param nb_iter_ecm : nombre d'itération pour le mélange du BG
		 * @param delta : paramètre de régularisation de l'algo de Mouloud
		 * @brief
		 * Constructeur
		 * 
		 */
		c_eyelids_segmentation ( 	unsigned int width,
									unsigned int height,
									unsigned int nb_samples,
									unsigned int nb_directions,
									unsigned int nb_iter_gem,
									unsigned int nb_iter_icm,
									unsigned int nb_iter_em,
									double delta,
									unsigned int closing,
									unsigned int opening,
									ostream * _err_stream = NULL
									);
		
		int default_setup(  unsigned int width,
							unsigned int height );
		
		/**@fn
		 * @param nb_samples : nombre d'échantillons radials
		 * @param nb_directions : nombre de directions
		 * @param width : largeur de l'image (width * r_ratio)
		 * @param height : hauteur de l'image (height * r_ratio )
		 * @param mean_bg : moyenne du bg
		 * @param sigma_bg : écart type du bg
		 * @param nb_iter_gem : nombre d'itération de l'algo de Mouloud
		 * @param nb_iter_icm : nombre d'itération de l'algo ICM
		 * @param nb_iter_ecm : nombre d'itération pour le mélange du BG
		 * @param delta : paramètre de régularisation de l'algo de Mouloud
		 * @brief
		 * Setup
		 * 
		 */	
		int setup( 	api_parameters & params, 
					unsigned int width,
					unsigned int height,
					ostream * _err_stream = NULL,
					const char * n_space = "eyelids", 
					const char * nb_directions_name = "nb_directions",
					const char * nb_samples_name = "nb_samples",
					const char * delta_name = "delta",
					const char * nb_iter_gem_name = "nb_iter_gem",
					const char * nb_iter_icm_name = "nb_iter_icm",
					const char * nb_iter_em_name = "nb_iter_em",
					const char * closing_name = "closing",
					const char * opening_name = "opening"
				);

		/**@fn
		 * @brief
		 * Setup
		 * 
		 */
		int setup ( void );
		
		/**@fn
		 * @param nb_samples : nombre d'échantillons radials
		 * @param nb_directions : nombre de directions
		 * @param width : largeur de l'image (width * r_ratio)
		 * @param height : hauteur de l'image (height * r_ratio )
		 * @param mean_bg : moyenne du bg
		 * @param sigma_bg : écart type du bg
		 * @param nb_iter_gem : nombre d'itération de l'algo de Mouloud
		 * @param nb_iter_icm : nombre d'itération de l'algo ICM
		 * @param nb_iter_ecm : nombre d'itération pour le mélange du BG
		 * @param delta : paramètre de régularisation de l'algo de Mouloud
		 * @brief
		 * Setup
		 * 
		 */	
		int setup ( unsigned int width,
					unsigned int height,
					unsigned int nb_samples,
					unsigned int nb_directions,
					unsigned int nb_iter_gem,
					unsigned int nb_iter_icm,
					unsigned int nb_iter_em,
					double delta,
					unsigned int closing,
					unsigned int opening,
					ostream * _err_stream = NULL
					);

		/**@fn
		 * @brief
		 * Destructeur.
		 * 
		 */
		~c_eyelids_segmentation();
		
		
		/**@fn
		 * @param polar_image : transformée polaire (normalisée) de l'iris
		 * @param polar_mask : masque
		 * @param x_pupil, y_pupil... : paramètres de la pupille
		 * @param x_iris, y_iris ... : paramètres de l'iris
		 * @param mean_pupil, sigma_pupil : approximation de la répartition des pixels de la pupille
		 * @brief
		 * Segmente la pupille
		 */
		int segment_eyelids ( 	const IplImage * polar_image,
								const IplImage * polar_mask,
								const IplImage * mask,
								unsigned int nb_samples_iris,
								const double & x_pupil,
								const double & y_pupil,
								const double * pupil_radii,
								const double * iris_radii,
								double mean_pupil,
								double sigma_pupil );
	
		inline unsigned int width() const
		{
			return _width;
		}
		
		inline unsigned int height() const
		{
			return _height;
		}
		
		inline const IplImage * polar_mask() const
		{
			return _polar_mask;
		}
		
		inline const IplImage * iris_mask() const
		{
			return _mask;
		}
		
		inline const double * iris_radii() const
		{
			return _radii;
		}
	protected:
		/**@fn
		 * @brief
		 * Lib. mémoire
		 * 
		 */
		void free();
		
		/**@fn
		 * @brief
		 * Ini. mémoire
		 * 
		 */
		void initialize();
		
		/**@fn
		 * @brief
		 * Alloc. mémoire
		 * 
		 */
		void alloc();
		
		/**@fn
		 * @brief
		 * Calcul du masque carthésien de l'iris 
		 * 
		 */
		void compute_masks(	const double * pupil_radii,
							const double * iris_radii,
							const double & x_pupil,
							const double & y_pupil,
							const IplImage * mask,
							const IplImage * polar_mask,
							unsigned int nb_samples_iris );
		
		Mouloud Mouloud_obj, Mouloud_obj_2;
		unsigned int  _width,
					   _height,
					   _nb_samples,
					   _nb_radii;
				
		unsigned int _nb_iter_gem,
					  _nb_iter_icm,
					  _nb_iter_em;
					  
		double _delta;			  
		IplImage * _mask;
		IplImage * _polar_mask;
		
		//Redimensionnement
		IplImage * _p_img,
				 * _p_mask;
		
		
		// two_gaussian_estimation_params est_bg;
		
		double * _radii;
		CvPoint * _points;
		IplConvKernel * structuring_element_1;
		IplConvKernel * structuring_element_2;
		unsigned int _closing, _opening;
		c_interpol_2d * inter;
		ostream * err_stream;
};


#endif
