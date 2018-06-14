/**@file c_pupil_thread.hpp
 * @brief
 * 
 * 
 */
#ifndef _C_PUPIL_THREAD_HPP_
#define _C_PUPIL_THREAD_HPP_
#include <unistd.h>
#include <iostream>
#include <opencv/cv.h>
#include <opencv/highgui.h>
#include "pupil_data.hpp"
#include "c_image_thread.hpp"
#include "c_preprocessing.hpp"
#include "c_pupil_tracking.hpp"
#include "c_pupil_segmentation.hpp"
using namespace std;


/**@class
 * @brief
 * Cette classe qui permet de gérer les étapes de preprocessing, de segmentation de la pupille et du calcul du premier score.
 * 
 */
class c_pupil_thread
{
	public:
		/**@fn
		 * @brief
		 * Constructeur
		 * 
		 */
		c_pupil_thread();
		
		/**@fn
		 * @brief
		 * Reset.
		 * 
		 */
		int setup ( void );
		
		int default_setup ( ostream * _err_stream = NULL );
		
		
		/**@fn
		 * @param argv : 
		 * - [0 => (n -1) ] : fichiers de paramètres
		 * @param argc : n
		 * @brief
		 * Setup.
		 * 
		 */
		int setup (	char ** argv,	
					unsigned int argc,
					ostream * _err_stream = NULL );	
		
		/**
		 * 
		 */
		int segment_pupil( 	const image_data & img_data );
						
		/**@fn
		 * @brief
		 * Destructeur.
		 */
		~c_pupil_thread();

		/**@fn
		 * @brief
		 * Reset des objets
		 * 
		 */
		void reset();


		/**@fn
		 * @return
		 * Données de segmentation de la pupille
		 * @warning 
		 * p_thread unsafe!
		 */
		inline const pupil_data & pupil_seg_data() const
		{
			return *p_data;
		}
		
		/**@fn
		 * @param error_str : flux d'erreur
		 * @brief
		 * Modifie le flux d'erreur.
		 * 
		 */
		inline void set_error_stream( ostream & error_str = cout )
		{
			err_stream = &error_str;
			preprocessing->set_error_stream( error_str );
			pupil_tracking->set_error_stream( error_str );
			pupil_segmentation->set_error_stream( error_str );
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
		
	
		//Obj.
		c_preprocessing * preprocessing;
		c_pupil_tracking * pupil_tracking;
		c_pupil_segmentation * pupil_segmentation;

		//Données
		pupil_data * p_data;
		

		

		//Flux d'erreur
		ostream * err_stream;
		
		int tracking;
};
#endif
