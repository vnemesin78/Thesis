/**@file tkalman_prediction.hpp
 * @author Valérian Némesin
 * @date 10/11/2011
 * @brief
 * Pairwise Kalman prediction
**/
#ifndef _TKALMAN_PREDICTION2_HPP_
	#define _TKALMAN_PREDICTION2_HPP_
	#include <gsl/gsl_matrix.h>
	#include <gsl/gsl_linalg.h>
	#include <gsl/gsl_blas.h>
	#include "gsl_triangle_matrix.hpp"

	/**@fn void tkalman_get_x_p(gsl_vector * x_p,
								const gsl_vector * _x_f,
								const gsl_vector * __y,
								const gsl_vector * _y,
								const gsl_matrix * f2_xx,
								const gsl_matrix * f2_xy,
								const gsl_matrix * q2_xy);
	 * @param[out] x_p : 
	 \f$ \hat{x}_{n|n - 1} \f$, predicting state expectation
	 * @param[in] _x_f : 
	 \f$ \hat{x}_{n - 1|n - 1} \f$, previous filtering state expectation
	 * @param[in] __y : 
	 \f$ y_{n-2} \f$, observation \f$(n-2)\f$
	 * @param[in] _y : 
	 \f$ y_{n-1} \f$, previous observation
	 * @param[in] f2_xx : 
	 \f$ F_2^{x,x} \f$, \f$ F^{x,x} - Q^{x,y} \: [Q^{y,y}]^{-1} \: F^{y,x} \f$
	 * @param[in] f2_xy : 
	 \f$ F_2^{x,y} \f$, \f$ F^{x,y} - Q^{x,y} \: [Q^{y,y}]^{-1} \: F^{y,y} \f$
	 * @param[in] q2_xy : 
	 \f$ Q_2^{x,y} \f$, \f$ Q^{x,y} \: [Q^{y,y}]^{-1} \f$
	 * @brief
	 This function computes the predicting state expectation:
	 - \f$ \hat{x}_{n|n - 1} = F_2^{x,x} \hat{x}_{n - 1|n - 1} + F_2^{x,y} y_{n-2} + Q_2^{x,y} y_{n-1}  \f$
	 */
	
	void tkalman_get_x_p(gsl_vector * x_p,
						 const gsl_vector * _x_f,
						 const gsl_vector * __y,
						 const gsl_vector * _y,
						 const gsl_matrix * f2_xx,
						 const gsl_matrix * f2_xy,
						 const gsl_matrix * q2_xy);
						 
	/**@fn void tkalman_get_sqrt_p_p_1(gsl_matrix * sqrt_p_p_1,
									   const gsl_matrix * sqrt_q_0_f,
									   const gsl_matrix * f2_xt,
									   const gsl_matrix * sqrt_q2_xx,
									   gsl_matrix * mat_xpt_x,
									   gsl_matrix * mat_xpt_x_view_00,
									   gsl_matrix * mat_xpt_x_view_10,
									   gsl_vector * vect_x)
	 * @param[out] sqrt_p_p_1 : 
	 \f$[P_{1|0}]^{\frac{1}{2}}\f$, square root of predicting state 1 covariance matrix
	 * @param[in] sqrt_q_0_f : 
	 \f$[Q_{0|0}]^{\frac{1}{2}}\f$, square root of filtering state 0 covariance matrix
	 * @param[in] f2_xt : 
	 \f$ F_2^{x,t} \f$, \f$ F^{x,t} - Q^{x,y} \: [Q^{y,y}]^{-1} \: F^{y,t} \f$
	 * @param[in] sqrt_q2_xx : 
	 \f$ [Q_2^{x,x}]^{\frac{1}{2}}\f$, square root of the reduced noise covariance.
	 * @param mat_xpt_x : 
	 \f$ M\f$,allocated \f$(n_x + n_t, n_x)\f$-matrix
	 * @param mat_xpt_x_view_00 : 
	 * view on the matrix \f$ M\f$, starting at \f$(0,0)\f$, ending at \f$(n_x - 1, n_x - 1)\f$
	 * @param mat_xpt_x_view_10 : 
	 * view on the matrix \f$ M\f$, starting at \f$(n_x,0)\f$, ending at \f$(n_x + n_t - 1, n_x - 1)\f$
	 * @param vect_x : allocated \f$(n_x)\f$-vector
	 * @brief
	 This function computes the square root of predicting state 1 covariance:
	 - 1. building the matrix
	 \f$
	 \begin{pmatrix}
		[Q_2^{x,x}]^{\frac{1}{2}} \\
		[Q_{0|0}]^{\frac{1}{2}} [F_2^{x,t}]'
	 \end{pmatrix}
	 \f$
	 - 2. QR decomposition, which can be written:
	 \f$
	 \begin{pmatrix}
		[P_{1|0}]^{\frac{1}{2}} \\
		0
	 \end{pmatrix}
	 \f$
	 * 
	 **/
	void tkalman_get_sqrt_p_p_1(gsl_matrix * sqrt_p_p_1,
								const gsl_matrix * sqrt_q_f_0,
								const gsl_matrix * f2_xt,
								const gsl_matrix * sqrt_q2_xx,
								gsl_matrix * mat_xpt_x,
								gsl_matrix * mat_xpt_x_view_00,
								gsl_matrix * mat_xpt_x_view_10,
								gsl_vector * vect_x);
								
	/**@fn void tkalman_get_sqrt_p_p(gsl_matrix * sqrt_p_p,
									 const gsl_matrix * _sqrt_p_f,
									 const gsl_matrix * f2_xx,
									 const gsl_matrix * sqrt_q2_xx,
									 gsl_matrix * mat_2xx,
									 gsl_matrix * mat_2xx_view_00,
									 gsl_matrix * mat_2xx_view_10,
									 gsl_vector * vect_x)
	 * @param[out] sqrt_p_p : 
	 \f$[P_{n|n-1}]^{\frac{1}{2}}\f$, square root of predicting state covariance matrix
	 * @param[in] _sqrt_p_f :
	 \f$[P_{n-1|n-1}]^{\frac{1}{2}}\f$, square root of previous filtering state covariance matrix
	 * @param[in] f2_xx : 
	 \f$ F_2^{x,x} \f$, \f$ F^{x,x} - Q^{x,y} \: [Q^{y,y}]^{-1} \: F^{y,x} \f$
	 * @param[in] sqrt_q2_xx : 
	 \f$ [Q_2^{x,x}]^{\frac{1}{2}}\f$, square root of the reduced noise covariance.
	 * @param mat_2xx : 
	 \f$ M\f$,allocated \f$(2n_x, n_x)\f$-matrix
	 * @param mat_2xx_view_00 : 
	 * view on the matrix \f$ M\f$, starting at \f$(0,0)\f$, ending at \f$(n_x - 1, n_x - 1)\f$
	 * @param mat_2xx_view_10 : 
	 * view on the matrix \f$ M\f$, starting at \f$(n_x,0)\f$, ending at \f$(2 n_x - 1, n_x - 1)\f$
	 * @param vect_x : allocated \f$(n_x)\f$-vector
	 * @brief
	 This function computes the square root of predicting state covariance:
	 - 1. building the matrix
	 \f$
	 \begin{pmatrix}
		[Q_2^{x,x}]^{\frac{1}{2}} \\
		[P_{n-1|n-1}]^{\frac{1}{2}} [F_2^{x,x}]'
	 \end{pmatrix}
	 \f$
	 - 2. QR decomposition, which can be written:
	 \f$
	 \begin{pmatrix}
		[P_{n|n-1}]^{\frac{1}{2}} \\
		0
	 \end{pmatrix}
	 \f$

	 */
	void tkalman_get_sqrt_p_p(gsl_matrix * sqrt_p_p,
							  const gsl_matrix * _sqrt_p_f,
							  const gsl_matrix * f2_xx,
							  const gsl_matrix * sqrt_q2_xx,
							  gsl_matrix * mat_2xx,
							  gsl_matrix * mat_2xx_view_00,
							  gsl_matrix * mat_2xx_view_10,
							  gsl_vector * vect_x);

	/**@fn void tkalman_do_prediction ( gsl_vector * x_p,
									    gsl_matrix * sqrt_p_p,
									    const gsl_vector * _x_f,
									    const gsl_matrix * _sqrt_p_f,
									    const gsl_vector * __y,
									    const gsl_vector * _y,
									    const gsl_matrix * f2_xx,
									    const gsl_matrix * f2_xy,
									    const gsl_matrix * sqrt_q2_xx,
									    const gsl_matrix * q2_xy,
									    gsl_matrix * mat_2xx,
									    gsl_matrix * mat_2xx_view_00,
									    gsl_matrix * mat_2xx_view_10,
									    gsl_vector * vect_x);
	 * @param[out] x_p : 
	 \f$ \hat{x}_{n|n - 1} \f$, predicting state expectation
	 * @param[out] sqrt_p_p : 
	 \f$[P_{n|n-1}]^{\frac{1}{2}}\f$, square root of predicting state covariance matrix
	 * @param[in] _x_f : 
	 \f$ \hat{x}_{n - 1|n - 1} \f$, previous filtering state expectation
	 * @param[in] _sqrt_p_f :
	 \f$[P_{n-1|n-1}]^{\frac{1}{2}}\f$, square root of previous filtering state covariance matrix
	 * @param[in] __y : 
	 \f$ y_{n-2} \f$, observation \f$(n-2)\f$
	 * @param[in] _y : 
	 \f$ y_{n-1} \f$, previous observation
	 * @param[in] f2_xx : 
	 \f$ F_2^{x,x} \f$, \f$ F^{x,x} - Q^{x,y} \: [Q^{y,y}]^{-1} \: F^{y,x} \f$
	 * @param[in] f2_xy : 
	 \f$ F_2^{x,y} \f$, \f$ F^{x,y} - Q^{x,y} \: [Q^{y,y}]^{-1} \: F^{y,y} \f$
	 * @param[in] sqrt_q2_xx : 
	 \f$ [Q_2^{x,x}]^{\frac{1}{2}}\f$, square root of the reduced noise covariance.
	 * @param[in] q2_xy : 
	 \f$ Q_2^{x,y} \f$, \f$ Q^{x,y} \: [Q^{y,y}]^{-1} \f$
	 * @param mat_2xx : 
	 \f$ M\f$,allocated \f$(2n_x, n_x)\f$-matrix
	 * @param mat_2xx_view_00 : 
	 * view on the matrix \f$ M\f$, starting at \f$(0,0)\f$, ending at \f$(n_x - 1, n_x - 1)\f$
	 * @param mat_2xx_view_10 : 
	 * view on the matrix \f$ M\f$, starting at \f$(n_x,0)\f$, ending at \f$(2 n_x - 1, n_x - 1)\f$
	 * @param vect_x : allocated \f$(n_x)\f$-vector
	 * @brief
	 This function performs the prediction step of Pairwise Kalman filter.
	 - 1. building the matrix
	 \f$
	 \begin{pmatrix}
		[Q_2^{x,x}]^{\frac{1}{2}} \\
		[P_{n-1|n-1}]^{\frac{1}{2}} [F_2^{x,x}]'
	 \end{pmatrix}
	 \f$
	 - 2. QR decomposition, which can be written:
	 \f$
	 \begin{pmatrix}
		[P_{n|n-1}]^{\frac{1}{2}} \\
		0
	 \end{pmatrix}
	 \f$
	 - 3. Computing predicting state expectation
		\f$ \hat{x}_{n|n - 1} = F_2^{x,x} \hat{x}_{n - 1|n - 1} + F_2^{x,y} y_{n-2} + Q_2^{x,y} y_{n-1}  \f$
	**/
	void tkalman_do_prediction ( gsl_vector * x_p,
							  gsl_matrix * sqrt_p_p,
							  const gsl_vector * _x_f,
							  const gsl_matrix * _sqrt_p_f,
							  const gsl_vector * __y,
							  const gsl_vector * _y,
							  const gsl_matrix * f2_xx,
							  const gsl_matrix * f2_xy,
							  const gsl_matrix * sqrt_q2_xx,
							  const gsl_matrix * q2_xy,
							  gsl_matrix * mat_2xx,
							  gsl_matrix * mat_2xx_view_00,
							  gsl_matrix * mat_2xx_view_10,
							  gsl_vector * vect_x 
							 );
							 
	/**@fn void tkalman_do_prediction_1 ( gsl_vector * x_p_1,
								   gsl_matrix * sqrt_p_p_1,
								   const gsl_vector * x_f_0,
								   const gsl_vector * y_f_m1,
								   const gsl_matrix * sqrt_q_f_0,
								   const gsl_vector * y_0,
								   const gsl_matrix * f2_xt,
								   const gsl_matrix * f2_xx,
								   const gsl_matrix * f2_xy,
								   const gsl_matrix * sqrt_q2_xx,
								   const gsl_matrix * q2_xy,
								   gsl_matrix * mat_xpt_x,
								   gsl_matrix * mat_xpt_x_view_00,
								   gsl_matrix * mat_xpt_x_view_10,
								   gsl_vector * vect_x
								  );
	 * @param[out] x_p_1 : 
	 \f$ \hat{x}_{1|0} \f$, predicting state 1 expectation
	 * @param[out] sqrt_p_p_1 : 
	 \f$[P_{1|0}]^{\frac{1}{2}}\f$, square root of predicting state 1 covariance matrix
	 * @param[in] x_f_0 : 
	 \f$ \hat{x}_{0|0} \f$, filtering state 0 expectation
	 * @param[in] y_f_m1 : 
	 \f$ \hat{y}_{-1|0} \f$, 
	 "filtering state" 0 expectation
	 * @param[in] sqrt_q_0_f : 
	 \f$ [Q_{0|0}]^{\frac{1}{2}} \f$, square root of initial filtering state covariance
	 * @param[in] y_0 : 
	 \f$y_0\f$, initial observation
	 * @param[in] f2_xt : 
	 \f$ F_2^{x,t} \f$, \f$ F^{x,t} - Q^{x,y} \: [Q^{y,y}]^{-1} \: F^{y,t} \f$
	 * @param[in] f2_xx : 
	 \f$ F_2^{x,x} \f$, \f$ F^{x,x} - Q^{x,y} \: [Q^{y,y}]^{-1} \: F^{y,x} \f$
	 * @param[in] f2_xy : 
	 \f$ F_2^{x,y} \f$, \f$ F^{x,y} - Q^{x,y} \: [Q^{y,y}]^{-1} \: F^{y,y} \f$
	 * @param[in] sqrt_q2_xx : 
	 \f$ [Q_2^{x,x}]^{\frac{1}{2}}\f$, square root of the reduced noise covariance.
	 * @param[in] q2_xy : 
	 \f$ Q_2^{x,y} \f$, \f$ Q^{x,y} \: [Q^{y,y}]^{-1} \f$
	 * @param mat_xpt_x : 
	 \f$ M\f$,allocated \f$(n_x + n_t, n_x)\f$-matrix
	 * @param mat_xpt_x_view_00 : 
	 * view on the matrix \f$ M\f$, starting at \f$(0,0)\f$, ending at \f$(n_x - 1, n_x - 1)\f$
	 * @param mat_xpt_x_view_10 : 
	 * view on the matrix \f$ M\f$, starting at \f$(n_x,0)\f$, ending at \f$(n_x + n_t - 1, n_x - 1)\f$
	 * @param vect_x : allocated \f$(n_x)\f$-vector
	 * @brief
	 This function performs the first prediction step of Pairwise Kalman filter.
	 - 1. building the matrix
	 \f$
	 \begin{pmatrix}
		[Q_2^{x,x}]^{\frac{1}{2}} \\
		[Q_{0|0}]^{\frac{1}{2}} [F_2^{x,t}]'
	 \end{pmatrix}
	 \f$
	 - 2. QR decomposition, which can be written:
	 \f$
	 \begin{pmatrix}
		[P_{1|0}]^{\frac{1}{2}} \\
		0
	 \end{pmatrix}
	 \f$
	 - 3. Computing predicting state expectation
		\f$ \hat{x}_{1|0} = F_2^{x,t} \hat{t}_{0|0} + Q_2^{x,y} y_{0}  \f$
	 */
	void tkalman_do_prediction_1 ( gsl_vector * x_p_1,
								   gsl_matrix * sqrt_p_p_1,
								   const gsl_vector * x_f_0,
								   const gsl_vector * y_f_m1,
								   const gsl_matrix * sqrt_q_f_0,
								   const gsl_vector * y_0,
								   const gsl_matrix * f2_xt,
								   const gsl_matrix * f2_xx,
								   const gsl_matrix * f2_xy,
								   const gsl_matrix * sqrt_q2_xx,
								   const gsl_matrix * q2_xy,
								   gsl_matrix * mat_xpt_x,
								   gsl_matrix * mat_xpt_x_view_00,
								   gsl_matrix * mat_xpt_x_view_10,
								   gsl_vector * vect_x
								  );
	/**@class tkalman_prediction
	 * @author 
     * Valérian Némesin
	 * @brief
	 * Pairwise Kalman prediction
	 */
	class tkalman_prediction
	{
		public:
			/**@fn tkalman_prediction :: tkalman_prediction(const gsl_matrix * f2_xt,
														    const gsl_matrix * sqrt_q2_xx,
														    const gsl_matrix * q2_xy)
			 * @param[in] f2_xt : 
			 \f$ F_2^{x,t} \f$, \f$ F^{x,t} - Q^{x,y} \: [Q^{y,y}]^{-1} \: F^{y,t} \f$
			 * @param[in] sqrt_q2_xx : 
			 \f$ [Q_2^{x,x}]^{\frac{1}{2}}\f$, square root of the reduced noise covariance.
			 * @param[in] q2_xy : 
			 \f$ Q_2^{x,y} \f$, \f$ Q^{x,y} \: [Q^{y,y}]^{-1} \f$
			 */
			tkalman_prediction(const gsl_matrix * f2_xt,
							   const gsl_matrix * sqrt_q2_xx,
							   const gsl_matrix * q2_xy);
		
			/**@fn int tkalman_prediction :: setup(const gsl_matrix * f2_x,
												   const gsl_matrix * sqrt_q2_xx,
												   const gsl_matrix * q2_xy)
			 * @param[in] f2_xt : 
			 \f$ F_2^{x,t} \f$, \f$ F^{x,t} - Q^{x,y} \: [Q^{y,y}]^{-1} \: F^{y,t} \f$
			 * @param[in] sqrt_q2_xx : 
			 \f$ [Q_2^{x,x}]^{\frac{1}{2}}\f$, square root of the reduced noise covariance.
			 * @param[in] q2_xy : 
			 \f$ Q_2^{x,y} \f$, \f$ Q^{x,y} \: [Q^{y,y}]^{-1} \f$
			 * @return
			 * - 0 OK
			 * - 1 problem
			 * @brief
			 * Setup
			**/
			virtual int setup(const gsl_matrix * f2_x,
							  const gsl_matrix * sqrt_q2_xx,
							  const gsl_matrix * q2_xy);
		
			/**@fn tkalman_prediction :: ~ tkalman_prediction()
			 */
			virtual ~tkalman_prediction();
			
			/**@fn bool tkalman_prediction :: operator!() const;
			 * @return
			 * - 0 OK
			 * - 1 memory error
			 * @brief
			 * This method checks object memory.
			 */
			virtual bool operator!() const;
			
			
			/**@fn inline void tkalman_prediction :: do_prediction ( gsl_vector * x_p,
																	 gsl_matrix * sqrt_p_p,
																	 const gsl_vector * _x_f,
																	 const gsl_matrix * _sqrt_p_f,
																	 const gsl_vector * __y,
																	 const gsl_vector * _y)
			 * @param[out] x_p : 
			 \f$ \hat{x}_{n|n - 1} \f$, predicting state expectation
			 * @param[out] sqrt_p_p : 
			 \f$[P_{n|n-1}]^{\frac{1}{2}}\f$, square root of predicting state covariance matrix
			 * @param[in] _x_f : 
			 \f$ \hat{x}_{n - 1|n - 1} \f$, previous filtering state expectation
			 * @param[in] _sqrt_p_f :
			 \f$[P_{n-1|n-1}]^{\frac{1}{2}}\f$, square root of previous filtering state covariance matrix
			 * @param[in] __y : 
			 \f$ y_{n-2} \f$, observation \f$(n-2)\f$
			 * @param[in] _y : 
			 \f$ y_{n-1} \f$, previous observation
			 * @brief
			 This function performs the prediction step of Pairwise Kalman filter.
			 - 1. building the matrix
			 \f$
			 \begin{pmatrix}
				[Q_2^{x,x}]^{\frac{1}{2}} \\
				[P_{n-1|n-1}]^{\frac{1}{2}} [F_2^{x,x}]'
			 \end{pmatrix}
			 \f$
			 - 2. QR decomposition, which can be written:
			 \f$
			 \begin{pmatrix}
				[P_{n|n-1}]^{\frac{1}{2}} \\
				0
			 \end{pmatrix}
			 \f$
			 - 3. Computing predicting state expectation
				\f$ \hat{x}_{n|n - 1} = F_2^{x,x} \hat{x}_{n - 1|n - 1} + F_2^{x,y} y_{n-2} + Q_2^{x,y} y_{n-1}  \f$
			 */
			inline void do_prediction ( gsl_vector * x_p,
										gsl_matrix * sqrt_p_p,
										const gsl_vector * _x_f,
										const gsl_matrix * _sqrt_p_f,
										const gsl_vector * __y,
										const gsl_vector * _y)
			{
				tkalman_do_prediction ( x_p,
										sqrt_p_p,
										_x_f,
										_sqrt_p_f,
										__y,
										_y,
										&f2_xx,
										&f2_xy,
										_sqrt_q2_xx,
										_q2_xy,
										&mat_2xx,
										&mat_2xx_view_00,
										&mat_2xx_view_10,
										vect_x 
									   );
			}
								 
			/**@fn inline void tkalman_prediction :: do_prediction_1 ( gsl_vector * x_p_1,
																	   gsl_matrix * sqrt_p_p_1,
																	   const gsl_vector * x_f_0,
																	   const gsl_vector * y_f_m1,
																	   const gsl_matrix * sqrt_q_f_0,
																	   const gsl_vector * y_0)
			 * @param[out] x_p_1 : 
			 \f$ \hat{x}_{1|0} \f$, predicting state 1 expectation
			 * @param[out] sqrt_p_p_1 : 
			 \f$[P_{1|0}]^{\frac{1}{2}}\f$, square root of predicting state 1 covariance matrix
			 * @param[in] x_f_0 : 
			 \f$ \hat{x}_{0|0} \f$, filtering state 0 expectation
			 * @param[in] y_f_m1 : 
			 \f$ \hat{y}_{-1|0} \f$, 
			 "filtering state" 0 expectation
			 * @brief
			 This function performs the first prediction step of Pairwise Kalman filter.
			 - 1. building the matrix
			 \f$
			 \begin{pmatrix}
				[Q_2^{x,x}]^{\frac{1}{2}} \\
				[Q_{0|0}]^{\frac{1}{2}} [F_2^{x,t}]'
			 \end{pmatrix}
			 \f$
			 - 2. QR decomposition, which can be written:
			 \f$
			 \begin{pmatrix}
				[P_{1|0}]^{\frac{1}{2}} \\
				0
			 \end{pmatrix}
			 \f$
			 - 3. Computing predicting state expectation
				\f$ \hat{x}_{1|0} = F_2^{x,t} \hat{t}_{0|0} + Q_2^{x,y} y_{0}  \f$
			**/
			inline void do_prediction_1 ( gsl_vector * x_p_1,
										  gsl_matrix * sqrt_p_p_1,
										  const gsl_vector * x_f_0,
										  const gsl_vector * y_f_m1,
										  const gsl_matrix * sqrt_q_f_0,
										  const gsl_vector * y_0)
			{
				tkalman_do_prediction_1 ( x_p_1,
										  sqrt_p_p_1,
										  x_f_0,
										  y_f_m1,
										  sqrt_q_f_0,
										  y_0,
										  _f2_x,
										  &f2_xx,
										  &f2_xy,
										  _sqrt_q2_xx,
										  _q2_xy,
										  mat_xpt_x,
										  &mat_xpt_x_view_00,
										  &mat_xpt_x_view_10,
										  vect_x
										 );
			}
			
			
		protected:
			/**@fn int tkalman_prediction :: set_params(const gsl_matrix * f2_x,
																const gsl_matrix * sqrt_q2_xx,
																const gsl_matrix * q2_xy)
			 * @param[in] f2_xt : 
			 \f$ F_2^{x,t} \f$, \f$ F^{x,t} - Q^{x,y} \: [Q^{y,y}]^{-1} \: F^{y,t} \f$
			 * @param[in] sqrt_q2_xx : 
			 \f$ [Q_2^{x,x}]^{\frac{1}{2}}\f$, square root of the reduced noise covariance.
			 * @param[in] q2_xy : 
			 \f$ Q_2^{x,y} \f$, \f$ Q^{x,y} \: [Q^{y,y}]^{-1} \f$
			 * @return
			 * - 0 OK
			 * - 1 problem
			 */
			int set_params(const gsl_matrix * f2_x,
						   const gsl_matrix * sqrt_q2_xx,
						   const gsl_matrix * q2_xy);
		
		
			/**@fn void tkalman_prediction :: free();
			 * @brief
			 * This function frees memory.
			 */
			void free();
		
			/**@fn int tkalman_prediction :: alloc();
			 * @return
			 * 0 si bon déroulement de l'op.
			 * @brief
			 * This function allocates object memory.
			 */
			int alloc();
		
			/**@fn void tkalman_prediction :: create_views();
			 * @brief
			 * This function creates matrix views.
			 */
			void create_views();
		
		
			/**@fn void tkalman_prediction :: initialize();
			 * @brief
			 * This function sets object attributes to 0.
			 */
			void initialize();
		//Données propres
			unsigned int _size_x;
			unsigned int _size_y;
			gsl_vector * vect_x; 
			gsl_matrix * mat_xpt_x;
			gsl_matrix mat_xpt_x_view_00;
			gsl_matrix mat_xpt_x_view_10;
			gsl_matrix mat_2xx;
			gsl_matrix mat_2xx_view_00;
			gsl_matrix mat_2xx_view_10;
			
			
		//Paramètres
			const gsl_matrix * _f2_x;
			gsl_matrix f2_xx;
			gsl_matrix f2_xy;
			const gsl_matrix * _sqrt_q2_xx;
			const gsl_matrix * _q2_xy;
	};
	
#endif
