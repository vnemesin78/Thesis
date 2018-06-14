/**@file c_image_thread.hpp
 * 
 */
#ifndef _C_IMAGE_THREAD_HPP_
#define _C_IMAGE_THREAD_HPP_
#include "c_mbgc_image_processing.hpp"
#include "c_preprocessing.hpp"
#include "image_data.hpp"
#include "lib_api.hpp"
#include <cstring>
#include <iostream>
#include <fstream>
#define DEFAULT_WIDTH 640
#define DEFAULT_HEIGHT 480
using namespace std;
/**@class
 * @brief
 * Gestion des prétraitements de l'image
 * 
 */
class c_image_thread
{
	public:
		/**@fn
		 * @brief
		 * Constructeur
		 * 
		 */
		c_image_thread( void );
		
		/**@fn
		 * @brief
		 * Setup. par défaut
		 * 
		 * 
		 */
		int default_setup ( ostream * _err_stream = NULL );
		
		/**@fn
		 * @brief
		 * Reset.
		 * 
		 */
		int setup ( void );
		
		/**@fn
		 * @brief
		 * Setup.
		 */
		int setup( char ** argv,
					unsigned int argc,	
					ostream * _err_stream = NULL );
		
		/**@fn
		 * @param image : image
		 * @brief
		 * Prétraitements
		 * 
		 */
		int process(	const IplImage * image, 
						const char * img_name);
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
		~c_image_thread();

		void reset() 
		{ 
			//_data->frame_id = 0;
		}

		void reset_id() 
		{ 
			_data->frame_id = 0;
		}

		inline const image_data & data() const
		{
			return *_data;
		}
		inline void set_score( const double & score ) const
		{
			_data->score = score;
		}
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
				

		image_data * _data;
		
		
		//Flux d'erreur
		ostream * err_stream;

		//Type ( 	1 <-> MBGC,
		//			0 <-> other )
		int type;
};
#endif
