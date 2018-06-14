/**@file iris_data.hpp
 * 
 */
#ifndef _IRIS_DATA_HPP_
	#define _IRIS_DATA_HPP_
	#include "pupil_data.hpp"
	
	/**@struct
	 * Paramètres pour l'initialisation de la structure de données de l'iris
	 * 
	 */
	struct iris_data_params
	{
		unsigned int 	img_width,
						img_height,
						polar_width,
						polar_height,
						nb_samples_iris,
						iris_code_width,
						iris_code_height,
						iris_width,
						iris_height;
	};



	/**@struct
	 * @brief
	 * Données de l'iris
	 * 
	 */
	struct iris_data
	{
		/**@fn
		 * @brief
		 * Constructeur
		 * 
		 */
		iris_data();
		
		/**@fn
		 * @param image_size : dimension de l'image
		 */
		iris_data( const iris_data_params & params ); 
		
		/**@fn
		 * @brief
		 * 
		 */
		iris_data( const iris_data & data ); 
		
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
		int setup ( const iris_data_params & params );
		
		/**@fn
		 * @brief
		 * Setup.
		 */
		int setup ( const iris_data & data );
		
		/**@fn
		 * @brief
		 * Op. de recopie.
		 * 
		 */
		iris_data & operator=( const iris_data & data );
		
		/**@fn
		 * @brief
		 * Destructeur
		 * 
		 */
		~iris_data();
		
		
		/**@fn
		 * @brief
		 * Comparaison des types
		 */
		bool operator==( const iris_data & data ) const;
		 
		/**@fn
		 * @brief
		 * Sauvegarde
		 * 
		 */
		int save ( const char * rep, unsigned int id ) const;
		
		int load ( const char * rep, unsigned int id );
		
		
		//Pupille
		pupil_data p_data;
		
		//Iris
		double x_iris,
				y_iris,
				a_iris,
				b_iris,
				theta_iris;
		
		double new_x_iris,
				new_y_iris,
				new_r_iris;
		
		//Transformée polaire normalisée
		IplImage * polar_image,
				 * polar_mask;
				
		IplImage * iris_image,
				 * iris_mask;
				
		unsigned int nb_samples;
		unsigned int nb_directions;
		unsigned int nb_samples_iris;
				 
		unsigned int 	iris_width,
						iris_height;
		double nrj_ratio;			 
		//Iris code
		unsigned int nb_samples_code;
		unsigned int nb_directions_code;
		IplImage * code,
				 * code_mask;
		
		//Données pour les paupières
		double a_upper,
			   c_upper,
			   theta_upper;
			   
		double a_lower,
			   c_lower,
			   theta_lower;		   
		
		
		
		bool seg_ok;
	};

	/**@fn
	 * @brief
	 * Alloc. pour le buffer
	 * 
	 */
	void * iris_data_alloc( const void * params );

	/**@fn
	 * @brief
	 * Copie pour le buffer
	 * 
	 */
	void iris_data_copy( 	void * data_tg,
							const void * data_src );

	/**@fn
	 * @brief
	 * Lib. mémoire pour le buffer
	 * 
	 */
	void iris_data_free( void * p );

	
	
	
	
#endif
