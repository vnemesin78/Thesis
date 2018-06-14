/**@file c_buffer.hpp
 * 
 */
#ifndef _C_BUFFER_HPP_
	#define _C_BUFFER_HPP_
	#include <sstream>
	#include <iostream>
	#include <exception> //Exceptions (pour ne pas faire planter le programme en cas de problèmes de mémoire ou autre)
	#include <stdexcept>
	#include <cstring>
	#include "c_buffer_data.hpp"
	using namespace std;

	
	/**@class c_buffer
	 */
	class c_buffer
	{
		public:
		
			/**@fn
			 * @brief
			 * Constructeur
			 * 
			 */
			c_buffer ( );
		
		
			/**@fn
			 * @param nb_obj : nombre d'objets
			 * @param alloc_function: fonction d'allocation
			 * @param copy_function : fonction de recopie
			 * @param free_function : fonction de lib_mémoire
			 * @param params : paramètres de construction
			 * @brief
			 * Constructeur
			 * 
			 */
			c_buffer ( unsigned int nb_obj,
						void * (*alloc_function) ( const void * params ),
						void (*copy_function) ( void * merde,
												const void * data ),
						void (*free_function) ( void * data ),
						const void * params,
						ostream * _err_stream = NULL );
			
			/**@fn
			 * @param nb_obj : nombre d'objets
			 * @param alloc_function: fonction d'allocation
			 * @param copy_function : fonction de recopie
			 * @param free_function : fonction de lib_mémoire
			 * @param params : paramètres de construction
			 * @brief
			 * Setup.
			 * 
			 */
			int setup ( unsigned int nb_obj,
						void * (*alloc_function) ( const void * params ),
						void (*copy_function) (	void * merde,
												const void * data ),
						void (*free_function) ( void * data ),
						const void * params,
					    ostream * _err_stream = NULL );
			
			
			/**@fn
			 * @brief
			 * Ajoute un objet.
			 */
			int add_object(	const c_buffer_data * data );
			
			int add_object( unsigned int id,
							double score,
							const void * data );
			
			
			/**@fn
			 * @param obj_id : id de l'objet
			 * @brief
			 * Supprime un objet de la liste
			 * 
			 */
			int delete_object( const c_buffer_data * data );
		
			/**@fn
			 * @param data : données lues
			 * @param score_id : id de l'objet
			 * @brief
			 * Lit un objet
			 **/
			int get_object (	c_buffer_data * data,
								unsigned int score_id );
		
			/**@fn
			 * @param data : données lues
			 * @param score_id : id de l'objet
			 * @brief
			 * Lit un objet
			 **/
			int get_object_with_id (	c_buffer_data * data,
										unsigned int obj_id );
			/**@fn
			 * @brief
			 * Renvoie si l'objet est inclu.
			 * 
			 **/
			int inclued ( const c_buffer_data * data );
		
			/**@fn
			 * @brief
			 * Renvoie le score min.
			 * 
			 **/
			double score_min();
		
			/**@fn
			 * @brief
			 * Destructeur
			 */
			~c_buffer();
		
			/**@fn
			 * @brief
			 * Reset des objets
			 * 
			 */
			void reset();		
		
			//Accesseur
			/**@fn
			 * @brief
			 * Renvoie le nombre d'objets max.
			 * 
			 */
			inline unsigned int nb_objects() const
			{
				return _nb_objects;
			}

		protected:
			/**@fn
			 * @brief
			 * free memory
			 */
			void free();
			
			/**@fn
			 * @brief
			 * Initialize object
			 * 
			 */
			void initialize();
			
			/**@fn
			 * @brief
			 * Efface l'obj. n°id.
			 * 
			 **/
			void erase_object( unsigned int id );
		
			/**@fn
			 * @brief
			 * Efface l'obj. n°id.
			 * 
			 **/
			void write_object(	unsigned int id,
								double score,
								const void * data,
								unsigned int obj_id );
		
			/**@fn
			 * @brief
			 * Efface l'obj. n°id.
			 * 
			 **/
			void read_object(	c_buffer_data * data,
								unsigned int id  );
		
		
			//Jetons
			int _score_editting,
				* _object_in_reading,
				* _object_in_writting;
		
		
			void start_object_writting ( unsigned int obj_id );
			void end_object_writting ( unsigned int obj_id );
			
			void start_object_reading ( unsigned int obj_id );
			void end_object_reading ( unsigned int obj_id );
			
			void start_score_editting ();
			void end_score_editting();
		
			//Nombre d'objets
			unsigned int _nb_objects;
			
			c_buffer_data * buffer_data;
			unsigned int * sorted_buffer_ids;

			ostream * err_stream;
	};
	
	
	
#endif
