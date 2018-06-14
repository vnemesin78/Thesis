/**@file moments.hpp
 * @brief
 * Ce fichier contient les fonctions qui permettent de calculer l'espérance et la variance d'une image. 
 * @todo
 * Moments d'ordre >
 * 
 **/
#ifndef _MOMENTS_HPP_
	#define _MOMENTS_HPP_
	#include <opencv/cv.h>
	#define GET_IMAGE_MEAN_TEMP(type) \
	template double get_image_mean ( const type * image_data,\
									 unsigned int width,\
									 unsigned int height,\
									 unsigned int width_step );
	
	#define GET_IMAGE_MEAN_TEMP_2(type) \
	template double get_image_mean ( const type * image_data,\
									 const unsigned char * mask_data,\
									 unsigned int width,\
									 unsigned int height,\
									 unsigned int width_step,\
									 unsigned int mask_width_step,\
									 unsigned char mask_value );
	
	#define GET_IMAGE_VARIANCE_TEMP(type) \
	template double get_image_variance ( const type * image_data,\
										 unsigned int width,\
										 unsigned int height,\
										 unsigned int width_step,\
										 double mean );
	
	#define GET_IMAGE_VARIANCE_TEMP_2(type) \
	template double get_image_variance ( const type * image_data,\
										 const unsigned char * mask_data,\
										 unsigned int width,\
										 unsigned int height,\
										 unsigned int width_step,\
										 unsigned int mask_width_step,\
										 unsigned char mask_value,\
										 double mean );
	
	
	
	
	
	
	/**@fn
	 * @param image_data : pixels de l'image
	 * @param width : largeur de l'image
	 * @param height : hauteur de l'image
	 * @param width_step : 
	 * Calcule la moyenne des pixels de l'image
	 * @return
	 * moyenne
	 * 
	 **/
	template<class type> double get_image_mean ( const type * image_data,
												 unsigned int width,
												 unsigned int height,
												 unsigned int width_step );
	
	/**@fn 
	 * @param image : image
	 * @brief
	 * Calcule la moyenne des pixels de l'image
	 * @return
	 * moyenne
	 **/
	double get_image_mean( const IplImage * image );



	/**@fn
	 * @param image_data : pixels de l'image
	 * @param width : largeur de l'image
	 * @param height : hauteur de l'image
	 * @param width_step : 
	 * Calcule la moyenne des pixels de l'image
	 * @return
	 * moyenne
	 * 
	 **/
	template<class type> double get_image_mean ( const type * image_data,
												 const unsigned char * mask_data,
												 unsigned int width,
												 unsigned int height,
												 unsigned int width_step,
												 unsigned int mask_width_step,
												 unsigned char mask_value );
	
	/**@fn 
	 * @param image : image
	 * @brief
	 * Calcule la moyenne des pixels de l'image
	 * @return
	 * moyenne
	 **/
	double get_image_mean( 	const IplImage * image,
							const IplImage * mask,
							unsigned char mask_value );

	/**@fn
	 * @param image_data : pixels de l'image
	 * @param width : largeur de l'image
	 * @param height : hauteur de l'image
	 * @param width_step :
	 * @param mean : moyenne de l'image 
	 * Calcule la variance des pixels de l'image
	 * @return
	 * Variance
	 * 
	 **/
	template<class type> double get_image_variance ( 	const type * image_data,
														unsigned int width,
														unsigned int height,
														unsigned int width_step,
														double mean );
	
	/**@fn
	 * @param image_data : pixels de l'image
	 * @param width : largeur de l'image
	 * @param height : hauteur de l'image
	 * @param width_step :
	 * @param mean : moyenne de l'image 
	 * Calcule la variance des pixels de l'image
	 * @return
	 * Variance
	 * 
	 **/
	template<class type> double get_image_variance ( 	const type * image_data,
														const unsigned char * mask_data,
														unsigned int width,
														unsigned int height,
														unsigned int width_step,
														unsigned int mask_width_step,
														unsigned char mask_value,
														double mean );
	
	/**@fn 
	 * @param image : image
	 * @param mean : moyenne de l'image 
	 * Calcule la variance des pixels de l'image
	 * @return
	 * Variance
	 **/
	double get_image_variance( 	const IplImage * image,
								double mean  );
	
	/**@fn 
	 * @param image : image
	 * @param mean : moyenne de l'image 
	 * Calcule la variance des pixels de l'image
	 * @return
	 * Variance
	 **/
	double get_image_variance( 	const IplImage * image,
								const IplImage * mask,
								unsigned char mask_value,
								double mean  );
	
	
	
	
#endif
