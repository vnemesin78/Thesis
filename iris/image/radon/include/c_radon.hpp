/**@file c_radon.hpp
 * @author Valerian
 * 
 **/
#ifndef _C_RADON_HPP_
	#define _C_RADON_HPP_
	#include <sstream>
	#include <iostream>
	#include <opencv/cv.h>
	#include <opencv/highgui.h>
	#include <exception> //Exceptions (pour ne pas faire planter le programme en cas de problèmes de mémoire ou autre)
	#include <stdexcept>
	#include <cstring>
	#include <cmath>
	using namespace std;
	
	
	/**@class
	 * @brief
	 * This class is required for Radon's transformation computation.
	 * 
	 * 
	 **/
	class c_radon
	{
		public:
			/**@fn
			 * @brief
			 * Constructor
			 * 
			 **/
			c_radon ( );
			
			/**@fn
			 * @param r_size_max : max radius
			 * @param nb_thetas_max : max number of thetas
			 * @brief
			 * Constructor
			 * 
			 **/
			c_radon ( unsigned int r_size_max,
					  unsigned int nb_thetas_max );
					  
	
			/**@fn
			 * @param r_size_max : max radius
			 * @param nb_thetas_max : max number of thetas
			 * @brief
			 * Setup
			 * 
			 **/
			int setup (	unsigned int r_size_max,
						unsigned int nb_theta_max );
		
		
			/**@fn
			 * @param img_data : image pixels
			 * @param img_width : image width
			 * @param img_height : image hieght
			 * @param img_width_step
			 * @param thetas : angles for the radon transform
			 * @param nb_thetas : number of angles
			 * @brief
			 * Compute the radon's transform of the image
			 *  
			 * 
			 **/
			template<class type> int compute_transform( 	const type * img_data, 
															unsigned int img_width,
															unsigned int img_height,
															unsigned int img_width_step,
															unsigned int nb_thetas );
			/**@fn
			 * @param image : image
			 * @param thetas : angles for the radon transform
			 * @param nb_thetas : number of angles
			 * @brief
			 * Compute the radon's transform of the image
			 * 
			 **/
			int compute_transform ( const IplImage * image,
									unsigned int nb_theta );
		

			/**@fn
			 * @brief
			 * Destructor
			 * 
			 **/
			 ~c_radon();
			 
			/**@fn
			 * @param r : radius
			 * @param theta : angle
			 * @brief
			 * Search the max of the radon transform.
			 * 
			 **/
			double search_max ( 	double & r,
									double & theta ) const;
							
			 
			//Accesseurs
			/**@fn
			 * @return
			 * Radon's transformation
			 * 
			 **/
			inline const double * radon_transform() const
			{
				return _radon_transform;
			}
			
			/**@fn
			 * @return
			 * width step
			 * 
			 **/
			inline unsigned int width_step() const
			{
				return _width_step;
			}
			
			/**@fn
			 * @return
			 * width step
			 * 
			 **/
			inline unsigned int r_size() const
			{
				return _r_size;
			}
			
			/**@fn
			 * @return
			 * width step
			 * 
			 **/
			inline unsigned int nb_thetas() const
			{
				return _nb_thetas;
			}
			
			
		protected:
		
			/**@fn
			 * @brief
			 * free memory
			 * 
			 **/
			void free();
		
			/**@fn
			 * @brief
			 * Ini. object.
			 * 
			 **/
			void initialize();
		
		
			unsigned int _r_size_max,
						 _nb_thetas_max,
						 _width_step;
			
			IplImage * radon_image;
			
			double * _radon_transform;
			unsigned int _r_size,
						 _nb_thetas;
		
			ostream * err_stream;
		
			//TMP
			double * x_table, //2 r_size_max
				   * y_table,
				   * x_cos_table,
				   * y_sin_table;
		
	};
	
#endif
