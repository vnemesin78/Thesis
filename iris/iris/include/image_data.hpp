/**@file image_data.hpp
 * @author Valérian Némesin
 * @brief
 * 
 */
#ifndef _IMAGE_DATA_HPP_
	#define _IMAGE_DATA_HPP_
	#include <opencv/cv.h>
	#include <opencv/highgui.h>
	#include "lib_api.hpp"
	#include "lib_image.hpp"
	#include <cstring>
	#include <iostream>
	#include <fstream>
	
	using namespace std;
	/**@struct
	 * @brief
	 * Données relatives à une image
	 * 
	 */
	struct image_data
	{
		/**@fn
		 * @brief
		 * Constructeur
		 * 
		 */
		image_data();
		
		/**@fn
		 * @param image_size : dimension de l'image
		 */
		image_data( const CvSize & image_size ); 
		
		/**@fn
		 * @brief
		 * 
		 */
		image_data( const image_data & data ); 
		
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
		 * Mise à jour des dimensions.
		 * 
		 */
		int set_dim( unsigned int w, unsigned int h);

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
		int setup ( const image_data & data );
		
		/**@fn
		 * @brief
		 * Op. de recopie.
		 * 
		 */
		image_data & operator=( const image_data & _p_data );
		
		/**@fn
		 * @brief
		 * Destructeur
		 * 
		 */
		~image_data();
		
		/**@fn
		 * @brief
		 * Comparaison des types
		 */
		bool operator==( const image_data & data ) const;
		
		/**@fn
		 * @brief
		 * Chargement
		 * 
		 */ 
		int load ( 	const char * rep,
					unsigned int id );
		
		
		int  save (	const char * rep,
					unsigned int id ) const;
		 
		 
		//Données
		IplImage * image;
		string name;		 
				 
		unsigned int width,
					   height;
		unsigned int frame_id;
		
		double score;
		
		//Flag de validité
		bool img_ok;
		
	};

	/**@fn
	 * @brief
	 * Alloc. pour le buffer
	 * 
	 */
	void * image_data_alloc( const void * params );

	/**@fn
	 * @brief
	 * Copie pour le buffer
	 * 
	 */
	void image_data_copy( 	void * data_tg,
							const void * data_src );

	/**@fn
	 * @brief
	 * Lib. mémoire pour le buffer
	 * 
	 */
	void image_data_free( void * p );
	
	
	
	
	
#endif
