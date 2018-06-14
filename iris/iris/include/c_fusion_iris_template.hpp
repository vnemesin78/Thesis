/**@file c_fusion_iris_template.hpp
 * 
 */
#ifndef _C_FUSION_IRIS_TEMPLATE_HPP_
	#define _C_FUSION_IRIS_TEMPLATE_HPP_
	#include <opencv/cv.h>
	#include <opencv/highgui.h>
	#include <cmath>
	#include "lib_image.hpp"
	#include "lib_api.hpp"
	#include "c_fusion.hpp"
	#include "c_focus_score.hpp"
	#include "distances.hpp"
	#include "iris_data.hpp"
	#include <iostream>
	#include <sstream>
	using namespace std;
	
	/**@class
	 * @brief
	 * Classe pour la gestion de la fusion d'iris code
	 * 
	 */
	class c_fusion_iris_template
	{
		public:
		
			/**@fn
			 * @brief
			 * Constructeur par défaut
			 */			
			 c_fusion_iris_template(  void	);
		
		
			/**@fn
			 * @param nb_images : nombre d'images fusionnées
			 * @param type : type de score 0 <-> lineaire 1 <-> exponentiel 2 <-> uniforme
			 * @param nb_directions : nombre de directions
			 * @param nb_samples : nombre de rayons
			 * @param d_theta : recadrage
			 * @param stream : flux d'erreur
			 * @brief
			 * Constructeur par défaut
			 */
			c_fusion_iris_template( 	unsigned int nb_images,
										int type,
										unsigned int nb_directions,
										unsigned int nb_samples,
										double d_theta,
										ostream * stream = NULL );
										
			/**@fn
			 * @param nb_images : nombre d'images fusionnées
			 * @param type : type de score 0 <-> lineaire 1 <-> exponentiel 2 <-> uniforme
			 * @param nb_directions : nombre de directions
			 * @param nb_samples : nombre de rayons
			 * @param d_theta : recadrage
			 * @param stream : flux d'erreur
			 * @brief
			 * Setup
			 * 
			 */
			int setup ( 	unsigned int nb_images,
							int type,
							unsigned int nb_directions,
							unsigned int nb_samples,
							double d_theta,
							ostream * stream = NULL );
		
			/**@fn
			 * @param params : paramètres
			 * @brief
			 * Setup
			 * 
			 */
			int setup ( 	api_parameters & params,
							const char * n_space = "fusion",
							const char * nb_images_n = "nb_images",
							const char * type_n = "type",
							const char * n_space_iris_code = "iris_code",
							const char * nb_dir_name = "nb_directions",
							const char * nb_samples_name = "nb_samples",	
							const char * d_theta_name = "d_theta",				
							ostream * stream = &cout );
		
			/**@fn
			 * @brief
			 * Fusionne à partir d'un répertoire
			 */
			int fusion ( const char * rep );
		
			/**@fn
			 * @brief
			 * Fusionne à partir d'un ensemble d'iris segmentées
			 */
			int fusion ( 	const iris_data * iris_data,
							unsigned int nb_iris );
		
		
			int save ( 	const char * rep,
						const char * code_n = "fuzzy_iris_code",
						const char * fragil_n = "fragility_map" );
		
			//~ int safe_save( 	const char * filename,
							//~ unsigned id,
							//~ const char * img_filename );
		
			/**@fn
			 * @brief
			 * Destructeur
			 * 
			 */
			~c_fusion_iris_template();
	
			/**@fn
			 * @param error_str : flux d'erreur
			 * @brief
			 * Modifie le flux d'erreur.
			 * 
			 */
			inline void set_error_stream( ostream & error_str = cout )
			{
				err_stream = &error_str;
				_fusion_obj->set_error_stream( error_str );;
				
			}
			
			int default_setup( );
		protected:
			
			void get_iris_code_and_fragility_map( IplImage * iris_code, IplImage * map );
			
			
			/**@fn
			 * @brief
			 * ajoute un iris dans la fusion.
			 */
			void add_iris( 	const IplImage * code,
							const IplImage * mask,
							double score );
			
			/**@fn
			 * @brief
			 * Réalise la fusion
			 */
			void fusion_data ( void );
			
			/**@fn
			 * @brief
			 * Ini.
			 * 
			 */
			void initialize();
			
			/**@fn
			 * @brief
			 * Lib. mémoire
			**/
			void free();
		
			/**@fn
			 * @brief Alloc. mémoire
			 */
			void alloc();
		
			//Flux d'erreurs
			ostream * err_stream;
		
		
			//Paramètres
			int _type;
			unsigned int 	_nb_directions,
							_nb_samples,
							_nb_images,
							_nb_images_max;
		
			//Iris code et leurs masques
			IplImage ** _iris_codes,
					 ** _iris_code_masks;
			double * _scores;		 
					 
					 
			//Objet de fusion
			c_fusion * _fusion_obj;

			//Différentes fusions
			IplImage	** _saved_codes,
						** _saved_fragility_maps;
						
			double d_theta;
	};



#endif
