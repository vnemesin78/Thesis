/**@file c_focus_score.hpp
 * @brief
 * Ce fichier contient la classe qui permet de calculer les scores de focus.
 **/
#ifndef _C_FOCUS_SCORE_HPP_
	#define _C_FOCUS_SCORE_HPP_
	#include <opencv/cv.h>
	#include <opencv/highgui.h>
	#include <cstring>
	#include "lib_api.hpp"
	#include "lib_image.hpp"
	#include <sstream>
	#include <iostream>
	
	#define C_FOCUS_SCORE_GET_SCORE(type)\
	template double c_focus_score :: get_score ( 	const type * image,\
														const unsigned char * mask,\
														unsigned int width,\
														unsigned int height,\
														unsigned int width_step,\
														unsigned int mask_width_step,\
														int x,\
														int y,\
														unsigned int w,\
														unsigned int h,\
														unsigned int new_w,\
														unsigned int new_h );\
	template double c_focus_score :: get_score ( 	const type * image,\
													unsigned int width,\
													unsigned int height,\
													unsigned int width_step);
														
	
	using namespace std;
	/**@class
	 * @brief
	 * Cette classe permet de calculer le score de focus d'une image
	 * 
	 **/
	class c_focus_score
	{
		public:
			/**@fn 
			 * 
			 * @brief
			 * Constructeur par défaut.
			 * 
			 **/
			inline c_focus_score( void )
			{
				initialize();
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
			
			/**@fn 
			 * @param width : largeur de l'image
			 * @param height : hauteur de l'image
			 * @param kernel : noyau du filtre
			 * @param c : énergie de l'image pour un score de 50%
			 */
			c_focus_score(	unsigned int width,
							unsigned int height,
							const gsl_matrix * kernel,
							double c ,
							ostream * _err_stream = NULL);
			
			int default_setup(	void );
				   
			/**@fn 
			 * @param width : largeur de l'image
			 * @param height : hauteur de l'image
			 * @param kernel : noyau du filtre
			 * @param c : énergie de l'image pour un score de 50%
			 */
			int setup (	unsigned int width,
						unsigned int height,
						const gsl_matrix * kernel,
						double c,
						ostream * _err_stream = NULL );
			
			/**@fn
			 * @param params : variables pré chargées
			 * @param width : largeur
			 * @param height : hauteur
			 * @param prefix : préfixe des noms de variables
			 * @param c_name : nom de c.
			 * @brief
			 * Setup
			 * 
			 * 
			 **/
			int setup(	api_parameters & params,
						ostream * _err_stream = NULL,
						const char * prefix = "score",
						const char * width_name = "width",
						const char * height_name = "height",
						const char * kernel_name = "kernel",
						const char * c_name = "c" );
			
			
			/**@fn 
			 * 
			 * @param img_data 
			 * @param img_mask
			 * @param width
			 * @param height
			 * @param width_step
			 * @return
			 * Score de l'image
			 */
			 
			 template<class type> double get_score ( 	const type * image,
														const unsigned char * mask,
														unsigned int width,
														unsigned int height,
														unsigned int width_step,
														unsigned int mask_width_step,
														int x,
														int y,
														unsigned int w,
														unsigned int h,
														unsigned int new_w,
														unsigned int new_h );
			 
			 template<class type> double  get_score ( 	const type * image,
														unsigned int width,
														unsigned int height,
														unsigned int width_step );
			 
			double get_score (	const IplImage * image,
								const IplImage * mask,
								const CvRect & roi,
								unsigned int new_w,
								unsigned int new_h );
								
			double get_score (	const IplImage * image );
								
								
			/**@fn
			 * @brief
			 * Destructeur
			 * 
			 **/
			~c_focus_score();
		
			inline unsigned int width() const
			{
				return _width;
			}
			
			inline unsigned int height() const
			{
				return _height;
			}
		
		protected:
			void initialize();
			void free();
			
			double _c;
			CvMat * kernel;
		
			//Message d'erreurs
			ostream * err_stream;
			
			c_interpol_2d interpol_2d;
			IplImage * tmp_img,
					 * tmp_mask,
					 * r_img,
					 * r_mask,
					 * f_img;
					 
			unsigned int 	_width,
							_height;
	};
	
	
	
#endif
