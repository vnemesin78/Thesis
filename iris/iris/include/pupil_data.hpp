/**@file pupil_data.hpp
 * 
 * 
 */
#ifndef _PUPIL_DATA_HPP_
	#define _PUPIL_DATA_HPP_
	#include <iostream>
	#include <opencv/cv.h>
	#include <opencv/highgui.h>
	#include "image_data.hpp"
	
	using namespace std;
		
		
	/**@struct pupil_data
	 * @brief
	 * Données de la pupille
	 */
	struct pupil_data
	{
		/**@fn
		 * @brief
		 * Constructeur
		 * 
		 */
		pupil_data();
		
		/**@fn
		 * @param image_size : dimension de l'image
		 */
		pupil_data( const CvSize & image_size ); 
		
		/**@fn
		 * @brief
		 * 
		 */
		pupil_data( const pupil_data & data ); 
		
		/**@fn
		 * @brief
		 * Ini.
		 * 
		 */
		void initialize();

		/**@fn
		 * @brief
		 * Lib. mémoire
		 * 
		 */
		void free();
		
		/**@fn
		 * @brief
		 * Allocation mémoire.
		 * 
		 */
		void alloc();
		
		/**@fn
		 * @brief
		 * Reset.
		 * 
		 */
		int setup();
		
		/**@fn
		 * @brief
		 * Setup.
		 * 
		 */
		int setup ( const CvSize & image_size );
		
		/**@fn
		 * @brief
		 * Setup.
		 */
		int setup ( const pupil_data & data );
		
		/**@fn
		 * @brief
		 * Op. de recopie.
		 * 
		 */
		pupil_data & operator=( const pupil_data & _p_data );
		
		/**@fn
		 * @brief
		 * Destructeur
		 * 
		 */
		~pupil_data();
		
		
		/**@fn
		 * @brief
		 * Comparaison des types
		 */
		bool operator==( const pupil_data & data ) const;
		/**@fn
		 * @brief
		 * Mise à jour des dimensions.
		 * 
		 */
		int set_dim( unsigned int w, unsigned int h);
		
		/**@fn
		 * @brief
		 * Sauvegarde
		 * 
		 */ 
		int save ( const char * rep, unsigned int id ) const;
		
		/**@fn
		 * @brief
		 * Chargement
		 * 
		 */ 
		int load ( const char * rep, unsigned int id ); 
		
		
		
		
		
		
		
		
		
		
		
		image_data _img_data; //width, height, score0, frame_id, image
		
		IplImage  * smoothed_image;
		
		double x_pupil,
				y_pupil,
				a_pupil,
				b_pupil,
				theta_pupil,
				pupil_threshold,
				score;
				
		bool seg_ok;
		CvRect roi;
	};

	/**@fn
	 * @brief
	 * Alloc. pour le buffer
	 * 
	 */
	void * pupil_data_alloc( const void * params );

	/**@fn
	 * @brief
	 * Copie pour le buffer
	 * 
	 */
	void pupil_data_copy( 	void * data_tg,
							const void * data_src );

	/**@fn
	 * @brief
	 * Lib. mémoire pour le buffer
	 * 
	 */
	void pupil_data_free( void * p );

	
	
#endif
