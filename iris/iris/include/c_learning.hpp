#ifndef _C_LEARNING_HPP_
	#define _C_LEARNING_HPP_
	#include <iostream>
	#include <fstream>
	#include "c_get_iris_template.hpp"
	using namespace std;
	
	#define NB_PATH_MAX 250
	#define NB_CLASS_MAX 300
	
	
	
	/**@struct
	 * @brief
	 * Cette structure gère les paramètres d'une classe
	 * 
	 */
	struct iris_class_data
	{
		
		/**@fn
		 * @param name: nom de la classe
		 * @brief
		 * Constructeur
		 */
		iris_class_data( 	const string & name,
							unsigned int nb_img_max = NB_PATH_MAX );
		
		/**@fn
		 * @brief
		 * Ajout d'une image dans la classe
		 * @return
		 * -1 -> si classe pleine
		 * 1 -> si l'image n'appartient pas à la classe
		 * 0 -> sinon
		 */
		int add_image( 	const string & img_path );
		
		const string * operator[] ( unsigned int id ) const
		{
			if ( id >= nb_images )
				return NULL;
			else
				return (paths + id);
		}
		
		
		/**@fn
		 * @brief
		 * Destructeur
		 */
		~iris_class_data();
		
		string name; // L_***** ...
		char type;// L / R
		
		unsigned int nb_images;
		unsigned int nb_images_max;
		string * paths;

	};
	
	/**@class
	 * @brief
	 * Classe pour l'apprentissage
	 * 
	 */
	class c_learning
	{
		public:
			/**@fn
			 * @brief
			 * Constructeur
			 * 
			 */
			c_learning ( 	);
		
		
		
			/**@fn
			 * @brief
			 * Constructeur
			 * 
			 */
			c_learning ( 	int argc,
							char ** argv,
							unsigned int nb_classes_max = NB_CLASS_MAX,
							unsigned int nb_img_max = NB_PATH_MAX );
			/**@fn
			 * @brief
			 * Aide pour le fichier
			 * 
			 */
			unsigned int help( 	ostream & out, 
									unsigned int id ) const;
		
			/**@fn
			 * @brief
			 * Setup
			 * 
			 */
			int setup( 	int argc,
						char ** argv,
						unsigned int nb_classes_max = NB_CLASS_MAX,
						unsigned int nb_img_max = NB_PATH_MAX );
		
			/**@fn
			 * @brief
			 * Segmentation
			 * 
			 */
			int segment( );
		
			/**@fn
			 * @brief
			 * Destructeur
			 */
			~c_learning();
		protected:
			void free();
			
			void initialize();
			
			//Données
			struct iris_class_data ** data;
			unsigned int _nb_classes;
			unsigned int 	_nb_classes_max,
							_nb_img_max;
		
			string _save_rep;
			c_get_iris_template obj;
			
			//Log
			ofstream * log;
			
	};
	
	
	
	
#endif
