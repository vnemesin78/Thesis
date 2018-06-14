/**@file c_fusion.hpp
 * 
 **/
#ifndef _C_FUSION_HPP_
	#define _C_FUSION_HPP_
	#include "lib_image.hpp"
	#include <opencv/cv.h>
	#include <iostream>
	using namespace std;
	
	#define C_FUSION_TEMPLATE_ADD_IMAGE(type)\
	template void c_fusion :: add_image ( 	const type * img_data,\
											unsigned int img_width_step,\
											const unsigned char * mask_data,\
											unsigned int mask_width_step,\
											double score,\
											int d_x );
	
	#define C_FUSION_TEMPLATE_COMPUTE(type)\
	template void c_fusion :: compute (	type * img_data,\
										unsigned int img_width_step,\
										unsigned char * mask_data,\
										int mask_width_step,\
										bool binary_mask,\
										double threshold );
	
	
	
	
	/**@class c_fusion
	 * @brief
	 * Cette classe est construite pour la fusion d'image et/ou d'iris code
	 * 
	 **/
	class c_fusion
	{
		public:
			/**@fn
			 * @brief
			 * Constructeur par défaut
			 * 
			 **/
			c_fusion ( void );
		
			/**@fn
			 * @param width : largeur de l'image
			 * @param height : hauteur de l'image
			 * @brief
			 * Constructeur
			 * 
			 */
			c_fusion ( unsigned int width,
						unsigned int height,
						ostream * _err_stream = NULL );
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
			 * @brief
			 * Reset l'image
			*/
			void reset( void );
		
			/**@fn
			 * @param width : largeur de l'image
			 * @param height : hauteur de l'image
			 * @return 
			 * -1 en cas de problèmes
			 * @brief
			 * Constructeur
			 **/
			int setup (	unsigned int width,
						unsigned int height,
						ostream * _err_stream = NULL );
			/**@fn
			 * @param image : image
			 * @param mask : masque de l'image
			 * @param score : score de l'image
			 * @param d_x : recalage de l'image
			 * 
			 * @brief
			 * Ajoute une image dans la fusion
			 * 
			 **/
			int add_image (	const IplImage * image,
							const IplImage * mask,
							double score,
							int d_x = 0 );
			
			template <class type> void add_image (	const type * img_data,
													unsigned int img_width_step,
													const unsigned char * mask_data,
													unsigned int mask_width_step,
													double score,
													int d_x );
			
			/**@fn
			 * @param image : image de sortie
			 * @param mask : masque de sortie
			 * @param binary_mask : masque binaire?
			 * @param threshold : seuillage 
			 * @brief
			 * 
			 * 
			 **/
			int compute (	IplImage * image,
							IplImage * mask,
							bool binary_mask = true,
							double threshold = 0  );

			
			template <class type> void compute (	type * img_data,
													unsigned int img_width_step,
													unsigned char * mask_data,
													int mask_width_step,
													bool binary_mask = true,
													double threshold = 0.5  );
			
			
		
			/**@fn
			 * @brief
			 * Destructeur
			 * 
			 **/
			~c_fusion( void );
		protected:
		
			/**@fn 
			 * @brief
			 * Lib. mémoire
			 **/
			void free();
		
			/**@fn 
			 * @brief
			 * Ini. mémoire
			 **/
			void initialize();
		
			//Fusion d'images
			IplImage * _fuzzy_image,
					 * _coef_image;
			double _sum_coef;
			unsigned int 	_width,
							_height;
			
			ostream * err_stream;
		
		
	};
	
	
#endif
