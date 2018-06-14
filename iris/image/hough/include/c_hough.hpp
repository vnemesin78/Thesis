/**@file c_hough.hpp
 * 
 **/
#ifndef _C_HOUGH_HPP_
	#define _C_HOUGH_HPP_
	#include <opencv/cv.h>
	using namespace std;
	
	/**@class
	 * @brief
	 * Cette classe permet de gérer la transformée de Hough pour un cercle.
	 * 
	 **/
	class c_hough
	{
		public:
			
			/**@fn
			 * @brief
			 * Constructor
			 * 
			 **/
			c_hough();
				
			/**@fn
			 * @param width : largeur de l'image
			 * @param height : hauteur de l'image
			 * @param r_max : rayon max.
			 * @brief
			 * Setup
			 * 
			 **/
			c_hough(	unsigned int width,
						unsigned int height,
						unsigned int r_max,
						unsigned int thickness = 1 ); // hough_space, width, height
				
				

			/**@fn
			 * @param width : largeur de l'image
			 * @param height : hauteur de l'image
			 * @param r_max : rayon max.
			 * @brief
			 * Setup
			 * 
			 **/
			int  setup (	unsigned int width,
							unsigned int height,
							unsigned int r_max,
							unsigned int thickness = 1 ); // hough_space, width, height, r_step, _width_step
				
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
											unsigned int width_step = 0,
											unsigned int r_min = 1,
											unsigned int r_max = 0 );
		
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
											unsigned int width,
											unsigned int height,
											unsigned int r_min = 1,
											unsigned int r_max = 0 );
											
			/**@fn
			 * @param x, y, r : circle parameters
			 * @param r_min, r_max, width, height : aera of searching.
			 * @brief
			 * Find the best ciricle.
			 * 
			 */								
			unsigned int get_best_circle (	int & x,
											int & y,
											unsigned int & r,
											unsigned int r_min,
											unsigned int r_max,
											unsigned int width,
											unsigned int height );
			
			
			/**@fn 
			 * @brief
			 * Destructor
			 * 
			 */
			~c_hough();
		
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
			 * @param x, y : coordinate of pixel
			 * @param r_min : min_radius.
			 * @param r_max : max_radius.
			 * @brief
			 * Add a pixel on the hough transform.
			 **/
			void add_circles ( 	unsigned int y, 
								unsigned int x, 
								unsigned int r_min, 
								unsigned int r_max );
			
			/**@fn
			 * 
			 **/
			void generate_templates ( ); 

			ostream * err_stream;
			
			
			unsigned int * _hough_space;
			unsigned int _width,
						 _height,
						 _r_max,
						 _r_step,
						 _width_step;
			unsigned int _thickness;
			int ** _circles_points;
			unsigned int * _nb_circles_points;
	};
	
	
	
	
	
	
#endif
