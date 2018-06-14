
#ifndef _C_RECOGNITION_HPP_
	#define _C_RECOGNITION_HPP_
	#include "c_matching.hpp"
	
	/**@class
	 * @brief
	 * Classe de gestion de la reconnaissance d'iris.
	 * 
	 */
	class c_recognition
	{
		public:
			/**@fn
			 * @brief
			 * Constructeur
			 * 
			 */
			c_recognition();
			
			/**@fn
			 * @param argc: nombre d'arg.
			 * @param argv: arg.
			 * @param id : nombre d'images
			 * @param nb_thresholds : nombre de seuils testés pour la génération des courbes roc...
			 * @brief
			 * Constructeur
			 * 
			 */
			c_recognition(  	int argc, 
								char ** argv,
								const char * iris_code_fn,
								const char * fragility_map_fn,
								unsigned int nb_fused_images,
								unsigned int nb_thresholds  );
		
			/**@fn
			 * @param argc: nombre d'arg.
			 * @param argv: arg.
			 * @param id : nombre d'images
			 * @param nb_thresholds : nombre de seuils testés pour la génération des courbes roc...
			 * @brief
			 * Setup
			 * 
			 */
			int setup(  	int argc, 
							char ** argv,
							const char * iris_code_fn,
							const char * fragility_map_fn,
							unsigned int nb_fused_images,
							unsigned int nb_thresholds  );
		
			
			/**@fn
			 * @brief
			 * Matching
			 * 
			 */
			int match( 	double fbr, 
						double alpha = 0);
			
			/**@fn
			 * @brief
			 * Sauvegarde
			 * 
			 */
			int save( const char * rep ) const;
		
			/**@fn
			 * @brief
			 * Destructeur
			 */
			~c_recognition();
		
			unsigned int n_video() const
			{
				return _nb_video;
			}
			
			unsigned int nb_videos() const
			{
				return videos->nb_classes();
			}	
			
		
			unsigned int help ( ostream & out, 	
								unsigned int id ) const;
			
			/**@fn
			 * @param error_str : flux d'erreur
			 * @brief
			 * Modifie le flux d'erreur.
			 * 
			 */
			inline void set_error_stream( ostream & error_str = cout )
			{
				videos->set_error_stream( error_str );
				images->set_error_stream( error_str );
			}
			
			bool use_alpha() const
			{
				return ( images->get_distance() ==  (distance::function_prototype_bis) distance :: fragile_bit_distance ||
						 images->get_distance() ==  (distance::function_prototype_bis) distance :: Hamming_FBD );
			}

			bool use_fbr() const
			{
				return ( images->get_distance() ==  (distance::function_prototype_bis) distance :: Hamming 				||
						 images->get_distance() ==  (distance::function_prototype_bis) distance :: fragile_bit_distance ||
						 images->get_distance() ==  (distance::function_prototype_bis) distance :: Hamming_FBD );
			}



			
			
		protected:
			/**
			 * 
			 */
			void initialize();
			
			/**@fn
			 * @brief
			 * Lib mémoires
			 * 
			 * 
			 */
			void free();
			
			void get_roc_data();
			
			void get_top_rank();
			
			void sort_distances();
		
			c_database * videos;
			c_matching * images;
			
			
			unsigned int _id;
			double _nb_thresholds;
			
			
			
			
			
			unsigned int 	* truth,
							* orders,
							* top_rank;
							
			double * distances,
					* v_p, 
					* f_p, 
					* v_n, 
					* f_n,
					* thresholds;
		
			double * i_distances,
					* g_distances;
					
			double err;
			unsigned int _nb_video;
			
			
	};
	
	
	
#endif
