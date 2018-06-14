
#ifndef _FITTING_ELLIPSE_HPP_
	#define _FITTING_ELLIPSE_HPP_
	#include <gsl/gsl_matrix.h>
	#include <gsl/gsl_eigen.h>	
	/**@class fitting_ellipse
	 * @brief
	 * This class lets fit an sample of points with an ellipse. 
	 * 
	 */
	class fitting_ellipse
	{
		public:
			/**@fn fitting_ellipse ();
			 * @brief Constructor of the class.
			 */
			fitting_ellipse ();
			
			/**@fn ~fitting_ellipse ();
			 * @brief
			 * Destructor of the class.
			 */
			~fitting_ellipse ();
			
			/**@fn  template<class type> void  fit( const type * data,
													unsigned int n);
			 * @param data : points (x1, y1 ...)
			 * @param n : number of points
			 * @brief
			 * This function fits the sample of points with an ellipse.
			 */
			template<class type> void  fit( const type * data,
											unsigned int n,
											unsigned int nb_iter,
											const unsigned char * mask);
			/**@fn inline double x_center() const
			 * @brief
			 * This function returns the absciss of the center of the ellipse.
			 **/	
			inline double x_center() const
			{
				return _x_center;
			}
			
			/**@fn inline double y_center() const
			 * @brief
			 * This function returns the ordinate of the center of the ellipse.
			 **/	
			inline double y_center() const
			{
				return _y_center;
			}
			
			/**@fn inline double a() const
			 * @brief
			 * This function returns the first axis size.
			 **/	
			inline double a() const
			{
				return _a;
			}
			
			/**@fn inline double b() const
			 * @brief
			 * This function returns the second axis size.
			 **/	
			inline double b() const
			{
				return _b;
			}
			
			/**@fn inline double theta() const
			 * @brief
			 * This function returns the angle between the absciss axis and the first axis of ellipse. 
			 **/	
			inline double theta() const
			{
				return _theta;
			}
			
			
			
			
		protected:
			template<class type> void  normalize_data ( const type * data,
														unsigned int n,
														const unsigned char * mask);
			template<class type> void  compute_moments( const type * data,
														unsigned int n,
														const unsigned char * mask);
			void compute_singular_values ();
			void unnormalize();
			void compute_ellipse_params();
			void free();
			void initialize();
			
			gsl_matrix * m_matrix,
					   * c_matrix;
			gsl_eigen_genv_workspace * workspace;
			double m_x,
				   m_y,
				   s_x,
				   s_y;
			double moments_4[5],
				   moments_3[4],
				   moments_2[3],
				   moments_1[2],
				   par[6];
			double _x_center,
				   _y_center,
				   _a,
				   _b,
				   _theta;
			gsl_vector_complex * alpha;
			gsl_vector * beta;
			gsl_matrix_complex * e_vectors;
			gsl_vector * vector;
	};
	
#endif
