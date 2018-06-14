
#ifndef _C_GET_IRIS_TEMPLATE_HPP_
#define _C_GET_IRIS_TEMPLATE_HPP_
	#include "c_iris_thread.hpp"
	/**@class
	 * @brief
	 * Génération d'un template d'iris
	 * 
	 */
	class c_get_iris_template
	{
		public:
			
			c_get_iris_template();
			~c_get_iris_template();
			/**@fn
			 * @brief
			 * Setup.
			 */
			int setup( char ** argv,
						unsigned int argc,
						ostream * stream = NULL );
			
			unsigned int help( ostream & out, unsigned int id ) const;
			
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
				
			}
			
			
			/**@fn
			 * @brief
			 * Setup à la main sans programmes
			 * 
			 **/
			int default_setup ( ostream * stream = NULL );
		
			/**@fn
			 * @brief
			 * Segmentation
			 */
			int segment( const IplImage * image, const char * name );
			
		
			inline const iris_data & data() const
			{
				return iris_thread->iris_seg_data();
			}
			
			int safe_save ( 	const char * template_file, 
								const char * output_file,
								const char * img_filename );
			
			void reset_id()
			{
				image_thread->reset_id();
			}
			
		protected:
			c_image_thread * image_thread;
			c_pupil_thread * pupil_thread;
			c_iris_thread * iris_thread;
		
	};


#endif
