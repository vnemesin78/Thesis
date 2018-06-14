/**@file
 * @brief
 * Fait de l'interpolation bilinéaire d'une image
 * 
 **/

#ifndef _C_INTERPOL_2D_HPP_
	#define _C_INTERPOL_2D_HPP_

	#define C_INTERPOL_2D_INTERPOLATE(type) \
	template int c_interpol_2d :: interpolate ( 	const type * image_data,\
													const IplImage * mask,\
													unsigned int image_width_step,\
													unsigned int nb_iter );
													
	#define C_INTERPOL_2D_GET_IMAGE(type) \
	template int c_interpol_2d :: get_image ( 	type * image_data,\
												unsigned int width,\
												unsigned int height,\
												unsigned int width_step ) const;
	#include <opencv/cv.h>
	#include <iostream>
	using namespace std;										
	class c_interpol_2d
	{
		public:
			c_interpol_2d ()
			{
				initialize();
			}
		
		
			/**@fn
			 * @param width : largeur de l'image
			 * @param height : hauteur de l'image
			 * @param factors : facteur d'interpolation
			 * L , D, U, R
			 * @brief
			 * constructeur
			 *
			 **/
			c_interpol_2d (	unsigned int width,
							unsigned int height,
							const float * factors );
		
			/**@fn
			 * @param width : largeur de l'image
			 * @param height : hauteur de l'image
			 * @param factors : facteur d'interpolation
			 * L , D, U, R
			 * @brief
			 * setup
			 **/
			int setup (	unsigned int width,
						unsigned int height,
						const float * factors  );
		
			template<class type> int interpolate ( 	const type * image_data,
													const IplImage * mask,
													unsigned int image_width_step,
													unsigned int nb_iter = (unsigned int) -1 );
		
			int interpolate( 	const IplImage * image,
								const IplImage * mask,
								unsigned int nb_iter = (unsigned int) -1 );
		
		
			/**@fn
			 * @brief
			 * Destructeur
			 **/
			~c_interpol_2d();
		
			/**@fn
			 * @brief
			 * Renvoie l'image
			 * 
			 **/
			inline const IplImage * image()
			{
				return _image;
			}
			
			/**@fn
			 * @brief
			 * Renvoie le masque
			 * 
			 **/
			inline const IplImage * mask()
			{
				return _curr_mask;
			}
		
			/**@fn
			 * @brief
			 * Calcule l'image en niveaux de gris
			 * 
			 **/
			int get_image( IplImage * image ) const;
		
			/**@fn
			 * @brief
			 * 
			 * 
			 **/
			template <class type> int get_image ( 	type * image_data,
													unsigned int width,
													unsigned int height,
													unsigned int width_step ) const;
		
		
		
		protected:
		
			/**@fn
			 * @brief
			 * Lib. mémoire
			 * 
			 **/
			void free();
		
			/**@fn
			 * @brief
			 * Ini.
			 * 
			 **/
			void initialize();
			
			/**@fn 
			 * @brief
			 * Interpolation
			 * 
			 */
			int interpolate(unsigned int width, unsigned int height );
			
			
			//
			float _factors[4];
			
			//Dimensions
			unsigned int _width, 
						 _height;
			
			//Images
			IplImage * _image,
					 * _previous_mask,
					 * _curr_mask;
					 
			ostream * err_stream;
			//Dilatation
			IplConvKernel * structuring_element_1;
	};
	
	
#endif
