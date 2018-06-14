
#ifndef _DISPLAY_FUNCTION_HPP_
	#define _DISPLAY_FUNCTION_HPP_
	#include <opencv/cv.h>
	#include "iris_data.hpp"
	
	#define INT_RGB( r, g, b) 65536 * (r) + 256 * (g) + (b)
	
	/**@class
	 * @brief
	 * Cette classe gère l'affichage des info de seg. de l'image
	 * 
	 **/
	class c_display_image_data
	{
		public:
			
			/**@fn
			 * @brief
			 * Constructeur
			 * 
			 **/	
			c_display_image_data();
			
		
		
			/**@fn
			 * @brief
			 * Constructeur
			 * 
			 **/
			c_display_image_data( 	unsigned int img_width,
									unsigned int img_height, 
									unsigned int data_width,
									unsigned int data_heigth,
									int font_size,
									int font_color );
			
			/**@fn
			 * @brief
			 * Setup
			 * 
			 **/						
			int setup ( );
			
			
			/**@fn
			 * @brief
			 * Setup
			 * 
			 **/						
			int setup ( 	unsigned int img_width,
							unsigned int img_height, 
							unsigned int data_width,
							unsigned int data_heigth,
							int font_size,
							int font_color );
			
			/**@fn
			 * @brief
			 *
			 **/
			int display( const image_data & data );
			
			/**@fn
			 * @brief
			 * Destructeur
			 * 
			 */
			~c_display_image_data();
			
			inline const IplImage * seg_image() const
			{
				return _seg_image;
			}
			
			inline IplImage * data_image()
			{
				return _data_image;
			}
			
			inline int data_offset() const
			{
				return y;
			}
			
			
		protected:
			/**@fn
			 * @brief
			 * Lib. mémoire
			 * 
			 **/
			void free();
			/**@fn
			 * @brief
			 * Ini mémoire
			 * 
			 **/
			void initialize();
		
			//Données
			IplImage * _seg_image,
					 * _data_image;
			
			//Paramètres
			int font_color,
				font_size;
				
			int y;
	};
	
	/**@class
	 * @brief
	 * Cette classe gère l'affichage des info de seg. de l'image
	 * 
	 **/
	class c_display_pupil_data
	{
		public:
			/**@fn
			 * @brief
			 * Constructeur
			 * 
			 **/	
			c_display_pupil_data();
			
		
		
			/**@fn
			 * @brief
			 * Constructeur
			 * 
			 **/
			c_display_pupil_data( 	unsigned int img_width,
									unsigned int img_height, 
									unsigned int data_width,
									unsigned int data_heigth,
									int font_size,
									int font_color,
									int thickness,
									int pupil_color );
			
			/**@fn
			 * @brief
			 * Setup
			 * 
			 **/						
			int setup ( );
			
			
			/**@fn
			 * @brief
			 * Setup
			 * 
			 **/						
			int setup ( 	unsigned int img_width,
							unsigned int img_height, 
							unsigned int data_width,
							unsigned int data_heigth,
							int font_size,
							int font_color,
							int thickness,
							int pupil_color );
			
			/**@fn
			 * @brief
			 *
			 **/
			int display( const pupil_data & data );
			
			/**@fn
			 * @brief
			 * Destructeur
			 * 
			 */
			~c_display_pupil_data();
			
			inline const IplImage * image() const
			{
				return d_image->seg_image();
			}
			
			inline const IplImage * pupil_image() const
			{
				return _pupil_img;
			}
			
			inline IplImage * data_image()
			{
				return _data_image;
			}
			
			inline int data_offset() const
			{
				return y;
			}
			
			
			
		protected:
			/**@fn
			 * @brief
			 * Lib. mémoire
			 * 
			 **/
			void free();
			/**@fn
			 * @brief
			 * Ini mémoire
			 * 
			 **/
			void initialize();
			
			
			int font_color,
				font_size,
				pupil_color,
				thickness;
			IplImage * _pupil_img;
			IplImage * _data_image;
			c_display_image_data * d_image;
			int y;
		
	};
	
	/**@class
	 * @brief
	 * Cette classe gère l'affichage des info de seg. de l'image
	 * 
	 **/
	class c_display_iris_data
	{
		public:
			/**@fn
			 * @brief
			 * Constructeur
			 * 
			 **/	
			c_display_iris_data();
			
		
		
			/**@fn
			 * @brief
			 * Constructeur
			 * 
			 **/
			c_display_iris_data( 	unsigned int img_width,
									unsigned int img_height, 
									unsigned int data_width,
									unsigned int data_heigth,
									int font_size,
									int font_color,
									int thickness,
									int pupil_color,
									int color_upper_eyelid,
									int color_lower_eyelid,
									int color_iris,
									int color_mask, //Nice try
									int iris_code_width,
									int iris_code_height,
									int polar_img_width,
									int polar_img_height );
			
			/**@fn
			 * @brief
			 * Setup
			 * 
			 **/						
			int setup ( );
			
			
			/**@fn
			 * @brief
			 * Setup
			 * 
			 **/						
			int setup ( 	unsigned int img_width,
							unsigned int img_height, 
							unsigned int data_width,
							unsigned int data_heigth,
							int font_size,
							int font_color,
							int thickness,
							int pupil_color,
							int color_upper_eyelid,
							int color_lower_eyelid,
							int color_iris,
							int color_mask,
							int iris_code_width,
							int iris_code_height,
							int polar_img_width,
							int polar_img_height );
			
			/**@fn
			 * @brief
			 *
			 **/
			int display( const iris_data & data );
			
			/**@fn
			 * @brief
			 * Destructeur
			 * 
			 */
			~c_display_iris_data();
			
			inline const IplImage * image() const
			{
				return d_pupil->image();
			}
			
			inline const IplImage * pupil_image() const
			{
				return d_pupil->pupil_image();
			}
			
			inline const IplImage * data_image() const
			{
				return _data_image;
			}
			
			inline int data_offset() const
			{
				return y;
			}
			
			inline const IplImage * iris_image() const
			{
				return _iris_img;
			}
			
			inline const IplImage * polar_image() const
			{
				return _polar_img;
			}
			
			inline const IplImage * iris_code() const
			{
				return _code_img;
			}	
			
			
		protected:
			/**@fn
			 * @brief
			 * Lib. mémoire
			 * 
			 **/
			void free();
			/**@fn
			 * @brief
			 * Ini mémoire
			 * 
			 **/
			void initialize();
		
		
			int font_color,
				font_size,
				pupil_color,
				thickness,
				color_upper_eyelid,
				color_lower_eyelid,
				color_iris,
				color_mask;
			
			IplImage * _iris_img,
					 * _polar_img,
					 * _code_img;
					 
			IplImage * _data_image;
			
			c_display_pupil_data * d_pupil;
			int y;
	};
	
	
	
	
	
	
	
	

								

	
	
	
	
	
	
#endif
