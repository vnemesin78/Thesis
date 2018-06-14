/**@file c_hough_ellipse.hpp
 * 
 **/
#ifndef _C_HOUGH_ELLIPSE_HPP_
	#define _C_HOUGH_ELLIPSE_HPP_
	#include <sstream>
	#include <iostream>
	#include <opencv/cv.h>
	#include <exception> //Exceptions (pour ne pas faire planter le programme en cas de problèmes de mémoire ou autre)
	#include <stdexcept>
	#include <cstring>
	#include <cmath>
	using namespace std;
	
	class c_hough_ellipse
	{
		public:
			/**@fn
			 * @brief
			 * Constructor
			 * 
			 **/
			c_hough_ellipse();
				
			/**@fn
			 * @param width : largeur de l'image
			 * @param height : hauteur de l'image
			 * @param r_max : rayon max.
			 * @brief
			 * Setup
			 * 
			 **/
			c_hough_ellipse(	unsigned int nb_x,
								unsigned int nb_y,
								unsigned int nb_thetas,
								unsigned int r_min,
								unsigned int r_max,
								unsigned int thickness = 1); // hough_space, width, height
	
			/**@fn
			 * @brief
			 * Reset
			 * 
			 */
			void setup( void );
	
			int setup (	unsigned int nb_x,
						unsigned int nb_y,
						unsigned int nb_thetas,
						unsigned int r_min,
						unsigned int r_max,
						unsigned int thickness = 1);
	
			/**@fn 
			 * @param img_data : image pixels
			 * @param width : image width ( <= _width)
			 * @param height : image height ( <= _height )
			 * @param width_step : image widthStep
			 * @param r_min : min radius ( >= 1 ) 
			 * @param r_max : max radius ( <= _r_max )
			 * @brief
			 * Compute the hough transform.
			 **/
			int compute_hough_transform( 	const unsigned char * img_data,
											unsigned int width,
											unsigned int height,
											unsigned int width_step,
											int x_start,
											int x_end,
											int y_start,
											int y_end,											
											unsigned int r_min,
											unsigned int r_max );
	
			/**@fn
			 * @param contour : contour
			 * @param nb_pts_contour : nb pts of contour
			 * @param width : image width ( <= _width)
			 * @param height : image height ( <= _height )
			 * @param r_min : min radius ( >= 1 ) 
			 * @param r_max : max radius ( <= _r_max )
			 * @brief
			 * Compute the hough transform.
			 */
			int compute_hough_transform( 	const unsigned int * contour,
											unsigned int nb_pts_contour,
											int x_start,
											int x_end,
											int y_start,
											int y_end,											
											unsigned int r_min,
											unsigned int r_max );
		
		
		
			/**@fn
			 * @param contour : contour
			 * @param nb_pts_contour : nb pts of contour
			 * @param width : image width ( <= _width)
			 * @param height : image height ( <= _height )
			 * @param r_min : min radius ( >= 1 ) 
			 * @param r_max : max radius ( <= _r_max )
			 * @brief
			 * Compute the hough transform.
			 */
			unsigned int search_best_ellipse( 	int & x,
												int & y,
												unsigned int & a,
												unsigned int & b,
												double & theta,
												int x_start,
												int x_end,
												int y_start,
												int y_end,											
												unsigned int r_min,
												unsigned int r_max );
		
		
		
		
		
			/**@fn
			 * @brief
			 * Destructeur.
			 */
			~c_hough_ellipse();
		
		protected:
			/**@fn
			 * @brief
			 * Lib. mémoire
			 * 
			 **/
			void free();
			
			/**@fn
			 * @brief
			 * Ini. objet.
			 * 
			 **/
			void initialize();
		
			/**@fn
			 * 
			 **/
			void generate_templates ( ); 
		
			void add_ellipses ( 	unsigned int y,
									unsigned int x,
									int x_start,
									int x_end,
									int y_start,
									int y_end,											
									unsigned int r_min,
									unsigned int r_max );									
		
		
		
			unsigned int * _hough_space; //width * height *  ( r_max )² * nb_thetha 
			unsigned int 	_nb_x,
							_nb_y,
							_nb_r,
							_nb_thetas;
							
			unsigned int 	_r_min,
							_r_max;
			
			unsigned int _thickness;
			ostream * err_stream;
							
			int ** _template_ellipses; //  4 PI ( r_max )³ x nb_thetas
			unsigned int * _nb_ellipses_points; //(r_max )² x nb_thetas
			
	};
	



#endif
