/**@file C_MBGC_image_processing.hpp
 * 
 **/
#ifndef _C_MBGC_IMAGE_PROCESSING_HPP_
	#define _C_MBGC_IMAGE_PROCESSING_HPP_
	#include <opencv/cv.h>
	#include <iostream>
	#include "lib_image.hpp"
	using namespace std;
	/**@class c_MBGC_image_processing
	 * 
	 **/
	class c_MBGC_image_processing
	{
		public:
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
			 * @param width : largeur
			 * @param height : hauteur
			 * 
			 * 
			 **/
			c_MBGC_image_processing ( 	unsigned int width = 640,
											unsigned int height = 480, 
											ostream * _err_stream = NULL);
			
			/**@fn
			 * @param width : largeur
			 * @param height : hauteur
			 * 
			 * 
			 **/
			int setup ( 	unsigned int width = 640,
							unsigned int height = 480, 
							ostream * _err_stream = NULL );
			/**@fn
			 * @return
			 * Image
			 * 
			 **/
			inline const IplImage * image() const
			{
				return _image;
			}
		
			int process( const IplImage * src_image );
		
			/**@fn
			 * @brief
			 * Destructeur
			 **/
			~c_MBGC_image_processing();
		
		protected:
			//Détection des objets à 127
			c_label * label_obj;
			
		
			IplImage * _image;
			unsigned int _width,
						 _height;
		
			ostream * err_stream;
			/**@fn
			 * @brief
			 * Lib. mémoire
			 **/
			void free();
	
			/**@fn
			 * @brief
			 * Ini.
			 * 
			 **/
			void initialize();
	
	};
	
	
	
#endif
