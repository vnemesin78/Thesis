
#ifndef _C_MATCHING_HPP_
	#define _C_MATCHING_HPP_
	#include <opencv/cv.h>
	#include <opencv/highgui.h>
	#include <dirent.h>
	#include "lib_api.hpp"
	#include "distances.hpp"
	#include "c_database.hpp"
	
	
	/**@class c_matching
	 * 
	 */
	class c_matching : public c_database
	{
		
		public:
			/**@fn
			 * @brief
			 * Constructeur par défaut
			 * 
			 */
			c_matching ( void );
			
			/**@brief
			 * @param rep : répertoire de la base de données
			 * @param dist_function : distance
			 * @param fusion_file : fichier de fusion ( fusion.png);
			 * @param fusion_mask : fichier de masque de fusion ( mask.png)
			 * @param stream : flux d'erreur
			 * Constructeur
			**/
			c_matching ( 	const char * rep_name, 
							const char * iris_code_file,
							const char * fragility_map_file,
							distance::function_prototype_bis dist,
							double d_theta,
							ostream * stream = NULL);
			
			/**@brief
			 * @param rep : répertoire de la base de données
			 * @param dist_function : distance
			 * @param fusion_file : fichier de fusion ( fusion.png);
			 * @param fusion_mask : fichier de masque de fusion ( mask.png)
			 * @param stream : flux d'erreur
			 * Constructeur
			**/
			virtual int setup ( const char * rep_name, 
								const char * iris_code_file,
								const char * fragility_map_file,
								distance::function_prototype_bis dist,
								double d_theta,
								ostream * stream = NULL);
			
			/**@brief
			 * @param rep : répertoire de la base de données
			 * @param dist_function : distance
			 * @param fusion_file : fichier de fusion ( fusion.png);
			 * @param fusion_mask : fichier de masque de fusion ( mask.png)
			 * @param stream : flux d'erreur
			 * Constructeur
			**/
			virtual int setup ( 	const char * rep_name, 
									const char * iris_code_file,
									const char * fragility_map_file,
									api_parameters & params,
									ostream * stream = NULL,
									const char * n_space = "matching",
									const char * distance_name = "distance",
									const char * d_theta_name = "d_theta" );			
			
			
			/**@fn
			 * @brief Matching
			 * @return
			 * Meilleure classe
			 */
			unsigned int matching (	double * distances,
									unsigned int * orders,
									unsigned int & truth,
									const IplImage * iris_code,
									const IplImage * fragility_map,
									const double & fbt,
									const char * class_name,
									const char * video_name = NULL );
									
			/**@fn
			 * @brief Matching
			 * @return
			 * Meilleure classe
			 */
			unsigned int matching_opt (	double * distances,
										unsigned int * orders,
										unsigned int & truth,
										const unsigned long long * _iris_codes_bis,
										const unsigned long long * _masks,
										const unsigned long long * _fragile_bits,
										const char * class_name,
										const char * video_name = NULL );				
									
									
									
			/**@fn
			 * @brief
			 * Destructeur
			 * 
			 */
			~c_matching();
			
			
			/**@fn
			 * @brief
			 * Modifie le alpha.
			 * 
			 */
			int set_alpha( double alpha );
			
			/**@fn
			 * @brief
			 * Renvoie la distance utilisée.
			 **/
			inline distance::function_prototype_bis get_distance() const
			{
				return dist;
				
			}
			

		protected:
		
			/**@fn void initialize();
			 * @brief
			 * INi
			 */
			void initialize();
			
			/**@fn
			 * @brief
			 * Lib. mémoire
			 * 
			 */
			void free();
		
			//Distance employée
			distance::function_prototype_bis dist;
			distance::function_prototype_64b dist_bis;

			
			double d_theta;
			void * dist_params;
			

	};
	
	
	
#endif
