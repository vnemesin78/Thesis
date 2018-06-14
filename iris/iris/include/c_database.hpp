/**@fn
 * @author Valérian Némesin
 * @brief
 * Gestion des bases de données
 */
#ifndef _C_DATABASE_HPP_
	#define _C_DATABASE_HPP_
	#include <opencv/cv.h>
	#include <opencv/highgui.h>
	#include <dirent.h>
	#include "lib_api.hpp"
	#include "iris_default.hpp"
	#include "lib_image.hpp"
	/**@class
	 * @brief
	 * Classe pour gérer une base de données d'iris codes
	 * 
	 */
	class c_database
	{
		public:
			/**@fn
			 * @brief
			 * Constructeur par défaut
			 * 
			 */
			c_database ( void );
			
			
			/**@fn
			 * @param rep : répertoire de la base de données
			 * @param iris_code_file : nom du fichier de l'iris code
			 * @param fusion_mask : nom du fichier de la carte de fragilité
			 * @param stream : flux d'erreurs ( NULL pour le désactiver)
			 * @brief
			 * Constructeur par défaut
			 * 
			 */			
			c_database ( 	const char * rep,
							const char * iris_code_file,
							const char * fragility_map_file,
							ostream * stream = NULL );
			
			/**@fn
			 * @param rep : répertoire de la base de données
			 * @param iris_code_file : nom du fichier de l'iris code
			 * @param fusion_mask : nom du fichier de la carte de fragilité
			 * @param stream : flux d'erreurs ( NULL pour le désactiver)
			 * @brief
			 * Setup
			 * 
			 */			
			virtual int setup ( const char * rep,
								const char * iris_code_file,
								const char * fragility_map_file,
								ostream * stream = NULL );
		
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
			 * Destructeur
			 * 
			 */
			~c_database();
			
			/**@fn
			 * @brief
			 * Renvoie l'iris code n°i
			 * 
			 */
			const IplImage * iris_code( unsigned int id ) const;
			
			/**@fn
			 * @brief
			 * Renvoie le masque d'occlusion de l'iris code n°i
			 * 
			 */
			const IplImage * fragility_map( unsigned int id ) const;
			
			unsigned long long * iris_code_bis( unsigned int id) const;
			
			
			unsigned long long * mask( unsigned int id) const;
			
			unsigned long long * fragile_bits( unsigned int id) const;	
			
						
			/**@fn
			 * @brief
			 * Renvoie le nom n°i
			 * 
			 */
			const char * name( unsigned int id ) const;
			
			/**@fn
			 * @brief
			 * Renvoie le nom de la classe n°i
			 * 
			 */
			const char * class_name( unsigned int id ) const;
			
			/**@fn
			 * @brief
			 * Renvoie le nombre de classes
			 * 
			 */
			inline unsigned int nb_classes() const
			{
				return _nb_classes;
			}
		
			void compute_fragile_bit_threshold( double fragile_bit_rate );
			
			double get_fragile_bit_threshold ( unsigned int id ) const;
			
			
		protected:
			
			
			/**@fn
			 * @brief
			 * Lib. mémoire
			 * 
			 */
			void free();
			
			/**@fn
			 * @brief
			 * Ini. mémoire
			 * 
			 */			
			void initialize();
		
			/**@fn
			 * @brief
			 * Alloc. mémoire
			 **/
			void alloc();
		
			/**@fn
			 * @brief
			 * Renvoie le nombre de classes
			 * 
			 */
			unsigned int get_nb_classes(	const char * rep ) const;
		
			/**@fn
			 * @brief
			 * Chargment de la base
			 * 
			 */
			int load (	const char * rep,
						const char * iris_code_file,
						const char * fragility_map_file );
		
			/**@fn
			 * @brief
			 * Calcule le taux de bits fragile pour chaque seuil.
			 * 
			 **/
			void compute_fragility_rate( void );
			
			double get_fbt( double fbr, unsigned id );
			
			
		
			//Nombre de classes
			unsigned int _nb_classes;
			//Noms des classes
			string * _class_names;
			//Noms des échantillons
			string * _names;
			
			//Cartes de fragilités
			IplImage ** _iris_codes,
					 ** _fragility_maps;
			unsigned long long ** _iris_codes_bis,
							   ** _masks,
							   ** _fragile_bits;
			unsigned int width_step;
			//Flux d'erreurs
			ostream * err_stream;
			
			//Palette pour les bits fragiles
			double ** _fragility_rate;
			
			double * _fragile_bit_thresholds;
			double _previous_fragile_bit_rate;
			
			c_histogram hist;
			
	};




#endif
