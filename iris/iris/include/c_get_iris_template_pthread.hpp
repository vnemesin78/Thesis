/**@file c_get_iris_template_pthread.hpp
 * @author Valérian Némesin
 */

#ifndef _C_GET_IRIS_TEMPLATE_PTHREAD_HPP_
#define _C_GET_IRIS_TEMPLATE_PTHREAD_HPP_
#define BUFFER_NAMESPACE "buffer"
#include <pthread.h>
#include "c_iris_thread.hpp"
#include "c_buffer.hpp"
#include "c_focus_score.hpp"
#include "display_functions.hpp"
#include <ctime>
/**@class c_get_iris_template_pthread
 * @brief
 * Gestion de la segmentation multi_thread de l'oeil
 * 
 * 
 */
class c_get_iris_template_pthread
{
	public:
		/**@fn
		 * @brief
		 * Constructeur
		 */
		c_get_iris_template_pthread ( void );
	
	
		/**@fn
		 * @param
		 * argv[0] <-> exe name
		 * argv[1] <-> data
		 * argv[2] <-> save
		 * argv[3 - n] <-> params
		 * @brief
		 * Setup
		 * 
		 */
		int setup ( unsigned int argc, char ** argv , const char * save_dir,
					ostream * stream = NULL, bool display = false);
	
		/**@fn
		 * @brief
		 * HElp!
		 * 
		 */
		unsigned int help ( ostream & out, unsigned int id ) const;
	
		/**@fn
		 * @brief
		 * Lancement de l'objet.
		 */
		int run ( 	const char * filename,
					int & end );
	
		/**@fn
		 * @brief
		 * Destructeur.
		 * 
		 */
		~c_get_iris_template_pthread();
		
		int test()
		{
			iris_thread->reset();
		}
		/**@fn
		 * @brief
		 * Save.
		 * 
		 */
		int save ( const char * rep_name );
	
			/**@fn
		 * @param error_str : flux d'erreur
		 * @brief
		 * Modifie le flux d'erreur.
		 * 
		 */
		inline void set_error_stream( ostream & error_str = cout )
		{
			image_thread->set_error_stream( error_str );
			pupil_thread->set_error_stream( error_str );
			iris_thread->set_error_stream( error_str );
			//~ focus_score_pupil->set_error_stream( error_str );
			//~ focus_score_iris->set_error_stream( error_str );
		}
		
		
	protected:
		/**@fn
		 * @brief Ini.
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
		 * Acquistion et prétraitement de l'image
		 * 
		 */
		void image_acquistion ( );
	
		/**@fn
		 * @brief
		 * Segmentation de la pupille
		 * 
		 */
		void pupil_segmentation ( );
	
		/**@fn
		 * @brief
		 * Segmentation de l'iris
		 * 
		 */
		void iris_segmentation ( );
	
		/**@fn
		 * @brief
		 * Affichage
		 * 
		 **/
		void display();
	
		//Pour l'acquisition d'image
		CvCapture * video;
	
		//Différents buffer
		c_buffer  * buffer_image, // Score = id
				  * buffer_pupil, 
				  * buffer_iris;
		
		//Données tmp pour les différents buffer
		c_buffer_data * b_data_image,
					  * b_data_pupil,
					  * b_data_iris;
		
		//Obj. de traitement
		c_image_thread * image_thread;
		c_pupil_thread * pupil_thread;
		c_iris_thread * iris_thread;
		
		//Calcul des scores
		c_focus_score 	* focus_score_pupil,
							* focus_score_iris;
		unsigned int 	r_width,
						r_height;
					
		//p_thread
		pthread_t tab_threads[4]; 
		
		//Variable d'arrêt
		int * p_end;
		bool 	img_end,
				pupil_end,
				iris_end;

		const char * _filename;
		
		ofstream * fps_file;
		unsigned int video_id;
		
		unsigned int nb_abb,
					 nb_frames,
					 nb_p_frames,
					 nb_i_frames;
		
		
		friend void * image_function ( void * params );
		friend void * pupil_function ( void * params );
		friend void * iris_function ( void * params );
		
		IplImage 	* d_image, 
					* d_pupil, 
					* d_iris;
		pthread_mutex_t mutex1, mutex2, mutex3;
		
		c_display_image_data * image_display_obj;
		c_display_pupil_data * pupil_display_obj;
		c_display_iris_data * iris_display_obj;
		
		
			
		bool display_on;
		
		
		
};
#endif
