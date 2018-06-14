/**@file distances.hpp
 * @author Valérian Némesin
 * @brief
 * This file contains the Hamming distance functions.
 *
 */
#ifndef _HAMMING_DISTANCE_HPP_
	#define _HAMMING_DISTANCE_HPP_
	#include <opencv/cv.h>
	#include <opencv/highgui.h>
	#include "lib_image.hpp"
	#define ACC 1

	namespace distance
	{
		//Prototype des fonctions de matching
		/**@typedef
		 * @param[out] d : distance entre les deux iris-codes
		 * @param[in] code_1 : iris code 1
		 * @param[in] fragility_map_1 : carte de fragilité 1
		 * @param[in] code_2 : iris code 2
		 * @param[in] fragility_map_2 : carte de fragilité 2
		 * @param[in] theta : angle de recadrage
		 * @param[in] params : paramètres de la distance
		 * @return
		 * - 0 si matching réussi
		 * - 1 si échec (0 pixel à comparer!)
		 * @brief
		 * Prototype pour la comparaison de deux iris codes
		 *
		 */
		typedef int (*function_prototype) ( double & d,
											const IplImage * code_1,
											const IplImage * fragility_map_1,
											const IplImage * code_2,
											const IplImage * fragility_map_2,
											int theta,
											const void * params );
		/**@typedef
		 * @param[out] d : distance entre les deux iris-codes
		 * @param[in] code_1_data : pixel de l'iris code 1
		 * @param[in] fragility_map_1_data : pixel de la carte de fragilité 1
		 * @param[in] code_2_data : pixel de l'iris code 2
		 * @param[in] fragility_map_2_data : pixel de la carte de fragilité 2
		 * @param[in] width : nombre de directions angulaires de l'iris code
		 * @param[in] height : 2 fois le nombre de rayons de l'iris code
		 * @param[in] code_1_width_step : taille d'une ligne de l'iris code 1 dans la mémoire
		 * @param[in] fragility_map_1_width_step : taille d'une ligne de la carte de fragilité 1 dans la mémoire 
		 * @param[in] code_2_width_step : taille d'une ligne de l'iris code 2 dans la mémoire
		 * @param[in] fragility_map_2_width_step : taille d'une ligne de la carte de fragilité 2 dans la mémoire
		 * @param[in] theta : angle de recadrage
		 * @param[in] params : paramètres de la distance
		 * @return
		 * - 0 si matching réussi
		 * - 1 si échec (0 pixel à comparer!)
		 * @brief
		 * Prototype pour la comparaison de deux iris codes
		 *
		 */
		typedef int (*function_prototype_bis) ( double & d,
												const unsigned char * code_1_data,
												const unsigned char * fragility_map_1_data,
												const unsigned char * code_2_data,
												const unsigned char * fragility_map_2_data,
												unsigned int width,
												unsigned int height,
												unsigned int code_1_width_step,
												unsigned int fragility_map_1_width_step,
												unsigned int code_2_width_step,
												unsigned int fragility_map_2_width_step,
												int theta,
												const void * params   );
		
		
		/**@typedef
		 * @param[out] d : distance entre les deux iris-codes
		 * @param[in] code_1_data : pixel de l'iris code 1
		 * @param[in] fragility_map_1_data : pixel de la carte de fragilité 1
		 * @param[in] code_2_data : pixel de l'iris code 2
		 * @param[in] fragility_map_2_data : pixel de la carte de fragilité 2
		 * @param[in] width : nombre de directions angulaires de l'iris code
		 * @param[in] height : 2 fois le nombre de rayons de l'iris code
		 * @param[in] code_1_width_step : taille d'une ligne de l'iris code 1 dans la mémoire
		 * @param[in] fragility_map_1_width_step : taille d'une ligne de la carte de fragilité 1 dans la mémoire 
		 * @param[in] code_2_width_step : taille d'une ligne de l'iris code 2 dans la mémoire
		 * @param[in] fragility_map_2_width_step : taille d'une ligne de la carte de fragilité 2 dans la mémoire
		 * @param[in] theta : angle de recadrage
		 * @param[in] params : paramètres de la distance
		 * @return
		 * - 0 si matching réussi
		 * - 1 si échec (0 pixel à comparer!)
		 * @brief
		 * Prototype pour la comparaison de deux iris codes
		 *
		 */				
		typedef int (*function_prototype_8b) ( 	double & d,
												const unsigned char * code_1_data,
												const unsigned char * mask_1_data,
												const unsigned char * fb_1_data,
												const unsigned char * code_2_data,
												const unsigned char * mask_2_data,
												const unsigned char * fb_2_data,
												unsigned int width,
												unsigned int height,
												unsigned int width_step,
												int theta,
												const void * params,
												unsigned char * _code_2_data,
												unsigned char * _mask_2_data,
												unsigned char * _fb_2_data );
		
		
		/**@typedef
		 * @param[out] d : distance entre les deux iris-codes
		 * @param[in] code_1_data : pixel de l'iris code 1
		 * @param[in] fragility_map_1_data : pixel de la carte de fragilité 1
		 * @param[in] code_2_data : pixel de l'iris code 2
		 * @param[in] fragility_map_2_data : pixel de la carte de fragilité 2
		 * @param[in] width : nombre de directions angulaires de l'iris code
		 * @param[in] height : 2 fois le nombre de rayons de l'iris code
		 * @param[in] code_1_width_step : taille d'une ligne de l'iris code 1 dans la mémoire
		 * @param[in] fragility_map_1_width_step : taille d'une ligne de la carte de fragilité 1 dans la mémoire 
		 * @param[in] code_2_width_step : taille d'une ligne de l'iris code 2 dans la mémoire
		 * @param[in] fragility_map_2_width_step : taille d'une ligne de la carte de fragilité 2 dans la mémoire
		 * @param[in] theta : angle de recadrage
		 * @param[in] params : paramètres de la distance
		 * @return
		 * - 0 si matching réussi
		 * - 1 si échec (0 pixel à comparer!)
		 * @brief
		 * Prototype pour la comparaison de deux iris codes
		 *
		 */				
		typedef int (*function_prototype_16b) ( double & d,
												const unsigned short int * code_1_data,
												const unsigned short int * mask_1_data,
												const unsigned short int * fb_1_data,
												const unsigned short int * code_2_data,
												const unsigned short int * mask_2_data,
												const unsigned short int * fb_2_data,
												unsigned int width,
												unsigned int height,
												unsigned int width_step,
												int theta,
												const void * params,
												unsigned short int * _code_2_data,
												unsigned short int * _mask_2_data,
												unsigned short int * _fb_2_data );		
		
		
		/**@typedef
		 * @param[out] d : distance entre les deux iris-codes
		 * @param[in] code_1_data : pixel de l'iris code 1
		 * @param[in] fragility_map_1_data : pixel de la carte de fragilité 1
		 * @param[in] code_2_data : pixel de l'iris code 2
		 * @param[in] fragility_map_2_data : pixel de la carte de fragilité 2
		 * @param[in] width : nombre de directions angulaires de l'iris code
		 * @param[in] height : 2 fois le nombre de rayons de l'iris code
		 * @param[in] code_1_width_step : taille d'une ligne de l'iris code 1 dans la mémoire
		 * @param[in] fragility_map_1_width_step : taille d'une ligne de la carte de fragilité 1 dans la mémoire 
		 * @param[in] code_2_width_step : taille d'une ligne de l'iris code 2 dans la mémoire
		 * @param[in] fragility_map_2_width_step : taille d'une ligne de la carte de fragilité 2 dans la mémoire
		 * @param[in] theta : angle de recadrage
		 * @param[in] params : paramètres de la distance
		 * @return
		 * - 0 si matching réussi
		 * - 1 si échec (0 pixel à comparer!)
		 * @brief
		 * Prototype pour la comparaison de deux iris codes
		 *
		 */				
		typedef int (*function_prototype_32b) ( double & d,
												const unsigned long int * code_1_data,
												const unsigned long int * mask_1_data,
												const unsigned long int * fb_1_data,
												const unsigned long int * code_2_data,
												const unsigned long int * mask_2_data,
												const unsigned long int * fb_2_data,
												unsigned int width,
												unsigned int height,
												unsigned int width_step,
												int theta,
												const void * params,
												unsigned long int * _code_2_data,
												unsigned long int * _mask_2_data,
												unsigned long int * _fb_2_data );							
										
		
		
		
		/**@typedef
		 * @param[out] d : distance entre les deux iris-codes
		 * @param[in] code_1_data : pixel de l'iris code 1
		 * @param[in] fragility_map_1_data : pixel de la carte de fragilité 1
		 * @param[in] code_2_data : pixel de l'iris code 2
		 * @param[in] fragility_map_2_data : pixel de la carte de fragilité 2
		 * @param[in] width : nombre de directions angulaires de l'iris code
		 * @param[in] height : 2 fois le nombre de rayons de l'iris code
		 * @param[in] code_1_width_step : taille d'une ligne de l'iris code 1 dans la mémoire
		 * @param[in] fragility_map_1_width_step : taille d'une ligne de la carte de fragilité 1 dans la mémoire 
		 * @param[in] code_2_width_step : taille d'une ligne de l'iris code 2 dans la mémoire
		 * @param[in] fragility_map_2_width_step : taille d'une ligne de la carte de fragilité 2 dans la mémoire
		 * @param[in] theta : angle de recadrage
		 * @param[in] params : paramètres de la distance
		 * @return
		 * - 0 si matching réussi
		 * - 1 si échec (0 pixel à comparer!)
		 * @brief
		 * Prototype pour la comparaison de deux iris codes
		 *
		 */				
		typedef int (*function_prototype_64b) ( double & d,
												const unsigned long long * code_1_data,
												const unsigned long long * mask_1_data,
												const unsigned long long * fb_1_data,
												const unsigned long long * code_2_data,
												const unsigned long long * mask_2_data,
												const unsigned long long * fb_2_data,
												unsigned int width,
												unsigned int height,
												unsigned int width_step,
												int theta,
												const void * params,
												unsigned long long * _code_2_data,
												unsigned long long * _mask_2_data,
												unsigned long long * _fb_2_data );							
		
		#define DIST_ROTATE(type)\
		template int distance::rotate( 	type * code_out,\
										type * mask_out,\
										type * fb_data_out,\
										const type * code_in,\
										const type * mask_in,\
										const type * fb_data_in,\
										unsigned int width,\
										unsigned int height,\
										unsigned int width_step,\
										int theta );
		
				
		/**@fn
		 * @brief
		 * Effectue la rotation rapide d'un iris code.
		 * 
		 **/
		template <class type> int rotate( 	type * code_out,
											type * mask_out,
											type * fb_data_out,
											const type * code_in,
											const type * mask_in,
											const type * fb_data_in,
											unsigned int width,
											unsigned int height,
											unsigned int width_step,
											int theta );
									
		#define DIST_CONVERT(type)\
		template int distance:: convert( 	type * code_out,\
											type * mask_out,\
											type * fb_data_out,\
											unsigned int & width_step,\
											const unsigned char * code_data,\
											const unsigned char * fragility_map_data,\
											unsigned int width,\
											unsigned int height,\
											unsigned int code_width_step,\
											unsigned int fragility_map_width_step,\
											unsigned int fm_threshold );\
		template int distance:: convert(	type * code_out,\
											type * mask_out,\
											type * fb_data_out,\
											unsigned int & width_step,\
											const IplImage * code,\
											const IplImage * fragility_map,\
											unsigned int fm_threshold );
		
		
								
		/**@fn
		 * @brief
		 * Convertit un iris code image en iris code binaire
		 * 
		 **/								
		template <class type> int convert( 	type * code_out,
											type * mask_out,
											type * fb_data_out,
											unsigned int & width_step,
											const unsigned char * code_data,
											const unsigned char * fragility_map_data,
											unsigned int width,
											unsigned int height,
											unsigned int code_width_step,
											unsigned int fragility_map_width_step,
											unsigned int fm_threshold );							
											
		/**@fn
		 * @brief
		 * Convertit un iris code image en iris code binaire
		 * 
		 **/											
		template <class type> int convert(	type * code_out,
											type * mask_out,
											type * fb_data_out,
											unsigned int & width_step,
											const IplImage * code,
											const IplImage * fragility_map,
											unsigned int fm_threshold );
			
			
		#define DIST_DEBUG(type)\
		template int distance::debug_convert( 	const type * code,\
												const type * mask,\
												const type * fb,\
												unsigned int width_step,\
												unsigned int width,\
												unsigned int height );
		
										
		/**@fn
		 * @brief
		 * Affiche l'iris code convertit
		 * 
		 **/								
		template <class type> int debug_convert( const type * code,
												 const type * mask,
												 const type * fb,
												 unsigned int width_step,
												 unsigned int width,
												 unsigned int height );	
												 
												 
		#define DIST_HAMMING(type)\
		template int distance::Hamming_opt ( 	double & d,\
												const type * code_1_data,\
												const type * mask_1_data,\
												const type * fb_data_1,\
												const type * code_2_data,\
												const type * mask_2_data,\
												const type * fb_data_2,\
												unsigned int width,\
												unsigned int height,\
												unsigned int width_step,\
												int theta,\
												const void * params,\
												type * _code_2_data,\
												type * _mask_2_data,\
												type * _fb_data_2 );										 							
		/**@fn
		 * @brief
		 * Calcul rapide de la distance de Hamming
		 * 
		 **/									
		template <class type> int  Hamming_opt ( 	double & d,
													const type * code_1_data,
													const type * mask_1_data,
													const type * fb_data_1,
													const type * code_2_data,
													const type * mask_2_data,
													const type * fb_data_2,
													unsigned int width,
													unsigned int height,
													unsigned int width_step,
													int theta,
													const void * params,
													type * _code_2_data,
													type * _mask_2_data,
													type * _fb_data_2 );		
													
													
		#define DIST_FBD(type)\
		template int distance::fragile_bit_distance_opt ( 	double & d,\
															const type * code_1_data,\
															const type * mask_1_data,\
															const type * fb_data_1,\
															const type * code_2_data,\
															const type * mask_2_data,\
															const type * fb_data_2,\
															unsigned int width,\
															unsigned int height,\
															unsigned int width_step,\
															int theta,\
															const void * params,\
															type * _code_2_data,\
															type * _mask_2_data,\
															type * _fb_data_2 );							
													
																	
		/**@fn
		 * @brief
		 * Calcul rapide de la distance des bits fragiles
		 * 
		 **/								
		template <class type> int  fragile_bit_distance_opt ( 	double & d,
																const type * code_1_data,
																const type * mask_1_data,
																const type * fb_data_1,
																const type * code_2_data,
																const type * mask_2_data,
																const type * fb_data_2,
																unsigned int width,
																unsigned int height,
																unsigned int width_step,
																int theta,
																const void * params,
																type * _code_2_data,
																type * _mask_2_data,
																type * _fb_data_2 );										
			
		#define DIST_H_FBD(type)\
		template int distance :: Hamming_FBD_opt ( 	double & d,\
													const type * code_1_data,\
													const type * mask_1_data,\
													const type * fb_data_1,\
													const type * code_2_data,\
													const type * mask_2_data,\
													const type * fb_data_2,\
													unsigned int width,\
													unsigned int height,\
													unsigned int width_step,\
													int theta,\
													const void * params,\
													type * _code_2_data,\
													type * _mask_2_data,\
													type * _fb_data_2 );	
																
		/**@fn
		 * @brief
		 * Calcul rapide de la distance des bits fragiles
		 * 
		 **/												
		template <class type> int  Hamming_FBD_opt ( 	double & d,
														const type * code_1_data,
														const type * mask_1_data,
														const type * fb_data_1,
														const type * code_2_data,
														const type * mask_2_data,
														const type * fb_data_2,
														unsigned int width,
														unsigned int height,
														unsigned int width_step,
														int theta,
														const void * params,
														type * _code_2_data,
														type * _mask_2_data,
														type * _fb_data_2 );										
																											
																
																		
		#define DIST_ALL(type)\
				DIST_ROTATE(type)\
				DIST_CONVERT(type)\
				DIST_DEBUG(type)\
				DIST_FBD(type)\
				DIST_HAMMING(type)\
				DIST_H_FBD(type)
		//Distance de Hamming
		/**@struct
		 * @var fragile_bit_threshold_1 : seuil de bit fragiles 1
		 * @var fragile_bit_threshold_2 : seuil de bit fragiles 2
		 * @brief
		 * Paramètre pour la distance de Hamming avec prise en compte des bits fragiles
		 **/
		struct Hamming_parameters
		{
			double fragile_bit_threshold_1;
			double fragile_bit_threshold_2;
		};
		
		/**@fn
		 * @param[out] d : (p 2 double) distance entre les deux iris-codes
		 * @param[in] code_1 : iris code 1
		 * @param[in] fragility_map_1 : carte de fragilité 1
		 * @param[in] code_2 : iris code 2
		 * @param[in] fragility_map_2 : carte de fragilité 2
		 * @param[in] theta : angle de recadrage
		 * @param[in] params : (Hamming_parameters*) Paramètre pour la distance de Hamming avec prise en compte des bits fragiles
		 * @return
		 * - 0 si matching réussi
		 * - 1 si échec (0 pixel à comparer!)
		 * @brief
		 * Calcule la distance de Hamming entre 2 iris codes tout en prenant en compte le masque d'occlusion des bits fragiles.
		 *
		 */
		int Hamming ( double & d,
					  const IplImage * code_1,
					  const IplImage * fragility_map_1,
					  const IplImage * code_2,
					  const IplImage * fragility_map_2,
					  int theta,
					  const void * params );
		
		/**@fn
		 * @param[out] d : distance entre les deux iris-codes
		 * @param[in] code_1_data : pixel de l'iris code 1
		 * @param[in] fragility_map_1_data : pixel de la carte de fragilité 1
		 * @param[in] code_2_data : pixel de l'iris code 2
		 * @param[in] fragility_map_2_data : pixel de la carte de fragilité 2
		 * @param[in] width : nombre de directions angulaires de l'iris code
		 * @param[in] height : 2 fois le nombre de rayons de l'iris code
		 * @param[in] code_1_width_step : taille d'une ligne de l'iris code 1 dans la mémoire
		 * @param[in] fragility_map_1_width_step : taille d'une ligne de la carte de fragilité 1 dans la mémoire 
		 * @param[in] code_2_width_step : taille d'une ligne de l'iris code 2 dans la mémoire
		 * @param[in] fragility_map_2_width_step : taille d'une ligne de la carte de fragilité 2 dans la mémoire
		 * @param[in] theta : angle de recadrage
		 * @param[in] params : paramètres de la distance
		 * @return
		 * - 0 si matching réussi
		 * - 1 si échec (0 pixel à comparer!)
		 * @brief
		 * Calcule la distance de Hamming entre 2 iris codes tout en prenant en compte le masque d'occlusion des bits fragiles.
		 *
		 */
		int Hamming ( 	double & d,
						const unsigned char * code_1_data,
						const unsigned char * fragility_map_1_data,
						const unsigned char * code_2_data,
						const unsigned char * fragility_map_2_data,
						unsigned int width,
						unsigned int height,
						unsigned int code_1_width_step,
						unsigned int fragility_map_1_width_step,
						unsigned int code_2_width_step,
						unsigned int fragility_map_2_width_step,
						int theta,
						const void * params );
		
		
		//Distance des bits fragiles
		/**@fn
		 * @param[out] d : (p 2 double) distance entre les deux iris-codes
		 * @param[in] code_1 : iris code 1
		 * @param[in] fragility_map_1 : carte de fragilité 1
		 * @param[in] code_2 : iris code 2
		 * @param[in] fragility_map_2 : carte de fragilité 2
		 * @param[in] theta : angle de recadrage
		 * @param[in] params : (Hamming_parameters*) Paramètre pour la distance de Hamming avec prise en compte des bits fragiles
		 * @return
		 * - 0 si matching réussi
		 * - 1 si échec (0 pixel à comparer!)
		 * @brief
		 * Calcule la distance de Hamming entre 2 masques de bits fragiles
		 *
		 */
		int fragile_bit_distance ( 	double & d,
									const IplImage * code_1,
									const IplImage * fragility_map_1,
									const IplImage * code_2,
									const IplImage * fragility_map_2,
									int theta,
									const void * params );
	
		/**@fn
		 * @param[out] d : distance entre les deux iris-codes
		 * @param[in] code_1_data : pixel de l'iris code 1
		 * @param[in] fragility_map_1_data : pixel de la carte de fragilité 1
		 * @param[in] code_2_data : pixel de l'iris code 2
		 * @param[in] fragility_map_2_data : pixel de la carte de fragilité 2
		 * @param[in] width : nombre de directions angulaires de l'iris code
		 * @param[in] height : 2 fois le nombre de rayons de l'iris code
		 * @param[in] code_1_width_step : taille d'une ligne de l'iris code 1 dans la mémoire
		 * @param[in] fragility_map_1_width_step : taille d'une ligne de la carte de fragilité 1 dans la mémoire 
		 * @param[in] code_2_width_step : taille d'une ligne de l'iris code 2 dans la mémoire
		 * @param[in] fragility_map_2_width_step : taille d'une ligne de la carte de fragilité 2 dans la mémoire
		 * @param[in] theta : angle de recadrage
		 * @param[in] params : paramètres de la distance
		 * @return
		 * - 0 si matching réussi
		 * - 1 si échec (0 pixel à comparer!)
		 * @brief
		 * Calcule la distance de Hamming entre 2 masques de bits fragiles
		 *
		 */
		int fragile_bit_distance ( 	double & d,
									const unsigned char * code_1_data,
									const unsigned char * fragility_map_1_data,
									const unsigned char * code_2_data,
									const unsigned char * fragility_map_2_data,
									unsigned int width,
									unsigned int height,
									unsigned int code_1_width_step,
									unsigned int fragility_map_1_width_step,
									unsigned int code_2_width_step,
									unsigned int fragility_map_2_width_step,
									int theta,
									const void * params );
		
		//Mix
		/**@struct
		 * @var p_Hamming : seuils pour les bits fragiles
		 * @var alpha : alpha * Hamming + (1-alpha) * FBD
		 * @brief
		 * Paramètre pour la distance de Hamming avec prise en compte des bits fragiles
		 **/
		struct Hamming_FBD_parameters
		{
			Hamming_parameters p_Hamming;
			double alpha;
		};
		
		
		/**@fn
		 * @param[out] d : distance entre les deux iris-codes
		 * @param[in] code_1_data : pixel de l'iris code 1
		 * @param[in] fragility_map_1_data : pixel de la carte de fragilité 1
		 * @param[in] code_2_data : pixel de l'iris code 2
		 * @param[in] fragility_map_2_data : pixel de la carte de fragilité 2
		 * @param[in] width : nombre de directions angulaires de l'iris code
		 * @param[in] height : 2 fois le nombre de rayons de l'iris code
		 * @param[in] code_1_width_step : taille d'une ligne de l'iris code 1 dans la mémoire
		 * @param[in] fragility_map_1_width_step : taille d'une ligne de la carte de fragilité 1 dans la mémoire 
		 * @param[in] code_2_width_step : taille d'une ligne de l'iris code 2 dans la mémoire
		 * @param[in] fragility_map_2_width_step : taille d'une ligne de la carte de fragilité 2 dans la mémoire
		 * @param[in] theta : angle de recadrage
		 * @param[in] params : (Hamming_FBD_parameters*) Paramètre pour la distance de Hamming avec prise en compte des bits fragiles
		 * @return
		 * - 0 si matching réussi
		 * - 1 si échec (0 pixel à comparer!)
		 * @brief
		 * Calcul la somme de la distance de Hamming + FBD.
		 *
		 */
		int Hamming_FBD( 	double & d,
							const unsigned char * code_1_data,
							const unsigned char * fragility_map_1_data,
							const unsigned char * code_2_data,
							const unsigned char * fragility_map_2_data,
							unsigned int width,
							unsigned int height,
							unsigned int code_1_width_step,
							unsigned int fragility_map_1_width_step,
							unsigned int code_2_width_step,
							unsigned int fragility_map_2_width_step,
							int theta,
							const void * params );
		
		/**@fn
		 * @param[out] d : (p 2 double) distance entre les deux iris-codes
		 * @param[in] code_1 : iris code 1
		 * @param[in] fragility_map_1 : carte de fragilité 1
		 * @param[in] code_2 : iris code 2
		 * @param[in] fragility_map_2 : carte de fragilité 2
		 * @param[in] theta : angle de recadrage
		 * @param[in] params : (Hamming_FBD_parameters*) Paramètre pour la distance de Hamming avec prise en compte des bits fragiles
		 * @return
		 * - 0 si matching réussi
		 * - 1 si échec (0 pixel à comparer!)
		 * @brief
		 * Calcul la somme de la distance de Hamming + FBD.
		 */
		int Hamming_FBD(  double & d,
						  const IplImage * code_1,
						  const IplImage * fragility_map_1,
						  const IplImage * code_2,
						  const IplImage * fragility_map_2,
						  int theta,
						  const void * params );
		
		
		
		
		
		//Espérance de la distance de Hamming
		
		/**@fn
		 * @param[out] d : (p 2 double) distance entre les deux iris-codes
		 * @param[in] code_1 : iris code 1
		 * @param[in] fragility_map_1 : carte de fragilité 1
		 * @param[in] code_2 : iris code 2
		 * @param[in] fragility_map_2 : carte de fragilité 2
		 * @param[in] theta : angle de recadrage
		 * @param[in] params : (Hamming_parameters*) Paramètre pour la distance de Hamming avec prise en compte des bits fragiles
		 * @return
		 * - 0 si matching réussi
		 * - 1 si échec (0 pixel à comparer!)
		 * @brief
		 * Calcule l'espérance de la distance de Hamming entre 2 iris codes tout en prenant en compte le masque d'occlusion des bits fragiles.
		 *
		 */
		int Hamming_expectation ( double & d,
								  const IplImage * code_1,
								  const IplImage * fragility_map_1,
								  const IplImage * code_2,
								  const IplImage * fragility_map_2,
								  int theta,
								  const void * params = NULL );
		
		/**@fn
		 * @param[out] d : distance entre les deux iris-codes
		 * @param[in] code_1_data : pixel de l'iris code 1
		 * @param[in] fragility_map_1_data : pixel de la carte de fragilité 1
		 * @param[in] code_2_data : pixel de l'iris code 2
		 * @param[in] fragility_map_2_data : pixel de la carte de fragilité 2
		 * @param[in] width : nombre de directions angulaires de l'iris code
		 * @param[in] height : 2 fois le nombre de rayons de l'iris code
		 * @param[in] code_1_width_step : taille d'une ligne de l'iris code 1 dans la mémoire
		 * @param[in] fragility_map_1_width_step : taille d'une ligne de la carte de fragilité 1 dans la mémoire 
		 * @param[in] code_2_width_step : taille d'une ligne de l'iris code 2 dans la mémoire
		 * @param[in] fragility_map_2_width_step : taille d'une ligne de la carte de fragilité 2 dans la mémoire
		 * @param[in] theta : angle de recadrage
		 * @param[in] params : paramètres de la distance
		 * @return
		 * - 0 si matching réussi
		 * - 1 si échec (0 pixel à comparer!)
		 * @brief
		 * Calcule l'espérance de la distance de Hamming entre 2 iris codes tout en prenant en compte le masque d'occlusion des bits fragiles.
		 *
		 */
		int Hamming_expectation ( 	double & d,
									const unsigned char * code_1_data,
									const unsigned char * fragility_map_1_data,
									const unsigned char * code_2_data,
									const unsigned char * fragility_map_2_data,
									unsigned int width,
									unsigned int height,
									unsigned int code_1_width_step,
									unsigned int fragility_map_1_width_step,
									unsigned int code_2_width_step,
									unsigned int fragility_map_2_width_step,
									int theta,
									const void * params = NULL );
		
		//Recadrage
		/**
		 * @param[out] d : distance entre les deux iris-codes
		 * @param[out] theta : angle de recadrage entre les deux iris-codes
		 * @param[in] code_1 : iris code 1
		 * @param[in] fragility_map_1 : carte de fragilité 1
		 * @param[in] code_2 : iris code 2
		 * @param[in] fragility_map_2 : carte de fragilité 2
		 * @param[in] theta_min : angle de recadrage min
		 * @param[in] theta_max : angle de recadrage max
		 * @param[in] dist : distance
		 * @param[in] params : paramètres de la distance
		 * @return
		 * - 0 si matching réussi
		 * - 1 si échec (0 pixel à comparer!)
		 * @brief
		 * Calcule l'angle de recadrage et sa distance entre 2 iris-codes
		 *
		 */
		int registering (	double & d,
							int & theta,
							const IplImage * code_1,
							const IplImage * fragility_map_1,
							const IplImage * code_2,
							const IplImage * fragility_map_2,
							const void * params,
							function_prototype_bis dist,
							int theta_min,
							int theta_max );
								
		/**
		 * @param[out] d : distance entre les deux iris-codes
		 * @param[out] theta : angle de recadrage entre les deux iris-codes
		 * @param[in] code_1_data : pixel de l'iris code 1
		 * @param[in] fragility_map_1_data : pixel de la carte de fragilité 1
		 * @param[in] code_2_data : pixel de l'iris code 2
		 * @param[in] fragility_map_2_data : pixel de la carte de fragilité 2
		 * @param[in] width : nombre de directions angulaires de l'iris code
		 * @param[in] height : 2 fois le nombre de rayons de l'iris code
		 * @param[in] code_1_width_step : taille d'une ligne de l'iris code 1 dans la mémoire
		 * @param[in] fragility_map_1_width_step : taille d'une ligne de la carte de fragilité 1 dans la mémoire 
		 * @param[in] code_2_width_step : taille d'une ligne de l'iris code 2 dans la mémoire
		 * @param[in] fragility_map_2_width_step : taille d'une ligne de la carte de fragilité 2 dans la mémoire
		 * @param[in] theta_min : angle de recadrage min
		 * @param[in] theta_max : angle de recadrage max
		 * @param[in] dist : distance
		 * @param[in] params : paramètres de la distance
		 * @return
		 * - 0 si matching réussi
		 * - 1 si échec (0 pixel à comparer!)
		 * @brief
		 * Calcule l'angle de recadrage et sa distance entre 2 iris-codes
		 *
		 */				
		int registering ( 	double & d,
							int & theta,
							const unsigned char * code_1_data,
							const unsigned char * fragility_map_1_data,
							const unsigned char * code_2_data,
							const unsigned char * fragility_map_2_data,
							unsigned int width,
							unsigned int height,
							unsigned int code_1_width_step,
							unsigned int fragility_map_1_width_step,
							unsigned int code_2_width_step,
							unsigned int fragility_map_2_width_step,
							const void * params,
							function_prototype_bis dist,
							int theta_min,
							int theta_max );							

		int registering_64b(	double & d,
								int theta,
								const unsigned long long * code_1_data,
								const unsigned long long * mask_1_data,
								const unsigned long long * fb_1_data,
								const unsigned long long * code_2_data,
								const unsigned long long * mask_2_data,
								const unsigned long long * fb_2_data,
								unsigned int width,
								unsigned int height,
								unsigned int width_step,
								int theta_min,
								int theta_max,
								function_prototype_64b dist,
								const void * params,
								unsigned long long * _code_2_data,
								unsigned long long * _mask_2_data,
								unsigned long long * _fb_2_data
								 );	
	};

#endif
