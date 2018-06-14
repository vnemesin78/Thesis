/**@file tkalman_smoothing.hpp
*@author Valérian Némesin
*@brief
Pairwise Kalman smoothing
**/
#ifndef _TKALMAN_SMOOTHING2_HPP_
	#define _TKALMAN_SMOOTHING2_HPP_
	#include <gsl/gsl_matrix.h>
	#include <gsl/gsl_linalg.h>
	#include <gsl/gsl_blas.h>
	#include "gsl_triangle_matrix.hpp"

	/**@fn void tkalman_get_x_s(	gsl_vector * x_s,
									const gsl_vector * x_f,
									const gsl_vector * x_p_,
									const gsl_vector * x_s_,
									const gsl_matrix * gain)
	 * @param[out] x_s : 
	 if n > 0
	 \f$ \hat{x}_{n|N} \f$, smoothing state expectation
	 else
	 \f$ \hat{t}_{0|N} \f$, smoothing state expectation
	 * @param[in] x_f : 
	 if n > 0
	 \f$ \hat{x}_{n|n} \f$, filtering state expectation
	 else
	 \f$ \hat{t}_{0|0} \f$, filtering state expectation
	 * @param[in] x_p_ : 
	 \f$ \hat{x}_{n + 1|n} \f$, following predicting state expectation
	 * @param[in] x_s_ : 
	 \f$ \hat{x}_{n + 1|N} \f$, following smoothing state expectation
	 @param[in] gain : 
	 \f$ K_{n|N} \f$, smoothing gain
	 * @brief
	 This function is intended to compute the smoothing state expectation:
	 \f$  \hat{x}_{n|N} = \hat{x}_{n|n} + K_{n|N} \: (\hat{x}_{n + 1|N} - \hat{x}_{n+1|n} \f$ 
	 */
	void tkalman_get_x_s(gsl_vector * x_s,
						 const gsl_vector * x_f,
						 const gsl_vector * x_p_,
						 const gsl_vector * x_s_,
						 const gsl_matrix * gain);
						 
	/**@fn void tkalman_get_smoothing_gain_0(gsl_matrix * s_gain,
											 const gsl_matrix * sqrt_q_f,
											 const gsl_matrix * sqrt_p_p_,
											 const gsl_matrix * f2_xt,
											 gsl_matrix * mat_xx,
											 gsl_matrix * mat_tx,
											 gsl_permutation * perm_x)
	 * @param[out] s_gain : 
	 * \f$ K_{0|N} \f$, smoothing gain
	 * @param[in] sqrt_q_f : 
	 \f$[Q_{0|0}]^{\frac{1}{2}}\f$, square root of filtering state 0 covariance matrix
	 * @param[in] sqrt_p_p_ : 
	 \f$[P_{1|0}]^{\frac{1}{2}}\f$, square root of predicting state 1 covariance matrix
	 * @param[in] f2_xt :
	 \f$ F_2^{x,t} \f$, \f$ F^{x,t} - Q^{t,x} \: [Q^{y,y}]^{-1} \: F^{y,t} \f$
	 * @param mat_xx : allocated \f$(n_x, n_x)\f$-matrix.
	 * @param mat_tx : allocated \f$(n_t, n_x)\f$-matrix.
	 * @param perm_x : allocated \f$n_x\f$-permutation
	 * @brief
	 * This function computes Pairwise Kalman smoothing gain:
	 * \f$K_{0|N} = Q_{0|0} [F_2^{x,t}]' [P_{1|0}]^{-1}\f$
	 */
	void tkalman_get_smoothing_gain_0(gsl_matrix * s_gain,
									  const gsl_matrix * sqrt_q_f,
									  const gsl_matrix * sqrt_p_p_,
									  const gsl_matrix * f2_xt,
									  gsl_matrix * mat_xx,
									  gsl_matrix * mat_tx,
									  gsl_permutation * perm_x);



	/**@fn void tkalman_get_smoothing_gain(gsl_matrix * s_gain,
										   const gsl_matrix * sqrt_p_f,
										   const gsl_matrix * sqrt_p_p_,
										   const gsl_matrix * f2_xx,
										   gsl_matrix * mat_xx,
										   gsl_permutation * perm_x)
	 * @param[out] s_gain : 
	 * \f$ K_{n|N} \f$, smoothing gain
	 * @param[in] sqrt_p_f : 
	 \f$[P_{n|n}]^{\frac{1}{2}}\f$, square root of the current filtering state covariance matrix
	 * @param[in] sqrt_p_p_ : 
	 \f$[P_{n+1|n}]^{\frac{1}{2}}\f$, square root of the following predicting state covariance matrix
	 * @param[in] f2_xx :
	 \f$ F_2^{x,x} \f$, \f$ F^{x,x} - Q^{x,y} \: [Q^{y,y}]^{-1} \: F^{y,x} \f$
	 * @param mat_xx : allocated \f$(n_x, n_x)\f$-matrix.
	 * @param perm_x : allocated \f$n_x\f$-permutation
	 * @brief
	 * This function computes Pairwise Kalman smoothing gain:
	 * \f$K_{n|N} = P_{n|n} [F_2^{x,x}]' [P_{n+1|n}]^{-1}\f$
	 */
	void tkalman_get_smoothing_gain(gsl_matrix * s_gain,
									const gsl_matrix * sqrt_p_f,
									const gsl_matrix * sqrt_p_p_,
									const gsl_matrix * f2_xx,
									gsl_matrix * mat_xx,
									gsl_permutation * perm_x);
									
	/**@fn void tkalman_get_sqrt_p_s_and_c_s(	gsl_matrix * sqrt_p_s,
												gsl_matrix * c_s,
												const gsl_matrix * sqrt_p_f,
												const gsl_matrix * sqrt_p_s_,
												const gsl_matrix * f2_xx,
												const gsl_matrix * sqrt_q2_xx,
												const gsl_matrix * s_gain,
												gsl_matrix * mat_3x2x,
												gsl_matrix * mat_3x2x_view_00,
												gsl_matrix * mat_3x2x_view_01,
												gsl_matrix * mat_3x2x_view_10,
												gsl_matrix * mat_3x2x_view_11,
												gsl_matrix * mat_3x2x_view_20,
												gsl_matrix * mat_3x2x_view_21,
												gsl_vector * vect_2x)
	 * @param[out] sqrt_p_s :
	 \f$ [P_{n|N}]^{\frac{1}{2}}\f$, square root of the current smoothing state covariance if n > 0
	 \f$ [Q_{0|N}]^{\frac{1}{2}}\f$ else
	 * @param[out] c_s : 
	 \f$ [P_{n + 1|N}]^{\frac{1}{2}}\ \; K_{n|N}^T \f$
	 * @param[in] sqrt_p_f : 
	 \f$[P_{n|n}]^{\frac{1}{2}}\f$, square root of the current filtering state covariance matrix if n > 0
	 \f$[Q_{0|0}]^{\frac{1}{2}}\f$, square root of the current filtering state covariance matrix else
	 * @param[in]  sqrt_p_s_ : 
	 \f$ [P_{n + 1|N}]^{\frac{1}{2}}  \f$, square root of the following smoothing state covariance
	 * @param[in] f2_xx :
	 \f$ F_2^{x,x} \f$, \f$ F^{x,x} - Q^{x,y} \: [Q^{y,y}]^{-1} \: F^{y,x} \f$ if n > 0
	 \f$ F_2^{x,t} \f$, \f$ F^{x,t} - Q^{x,y} \: [Q^{y,y}]^{-1} \: F^{y,t} \f$ else
	 * @param[in] sqrt_q2_xx : 
	 \f$ [Q_2^{x,x}]^{\frac{1}{2}}\f$, square root of the reduced noise covariance.
	 * @param[in] s_gain : 
	 * \f$ K_{n|N} \f$, smoothing gain 
	 * @param mat_3x2x : 
	 - \f$M\f$, allocated \f$(3 n_x, 2 n_x)\f$-matrix if n > 0
	 - \f$M\f$, allocated \f$(2 n_x + n_t, n_x + n_t)\f$-matrix if n = 0
	 * @param mat_3x2x_view_00 : 
	 * matrix view on \f$M\f$, starting at \f$(0,0)\f$, ending at \f$(n_x - 1, n_x - 1)\f$
	 * @param mat_3x2x_view_01 : 
	 * - matrix view on \f$M\f$, starting at \f$(0,n_x)\f$, ending at \f$(n_x - 1, 2n_x - 1)\f$ if n > 0
	 * - matrix view on \f$M\f$, starting at \f$(0,n_x)\f$, ending at \f$(n_x - 1, n_x + n_t - 1)\f$ else
	 * @param mat_3x2x_view_10 : 
	 * - matrix view on \f$M\f$, starting at \f$(n_x,0)\f$, ending at \f$(2 n_x - 1, n_x - 1)\f$ if n > 0
	 * - matrix view on \f$M\f$, starting at \f$(n_x,0)\f$, ending at \f$(n_x + n_t - 1, n_x - 1)\f$ else
	 * @param mat_3x2x_view_11 : 
	 * - matrix view on \f$M\f$, starting at \f$(n_x,n_x)\f$, ending at \f$(2 n_x - 1, 2 n_x - 1)\f$ if n > 0
	 * - matrix view on \f$M\f$, starting at \f$(n_x,n_x)\f$, ending at \f$(n_x + n_t - 1, n_x + n_t - 1)\f$ else
	 * @param mat_3x2x_view_20 :
	 * - matrix view on \f$M\f$, starting at \f$(2 n_x,0)\f$, ending at \f$(3 n_x - 1, n_x - 1)\f$ if n > 0
	 * - matrix view on \f$M\f$, starting at \f$(n_x + n_t,0)\f$, ending at \f$(2n_x + n_t - 1, n_x - 1)\f$ else
	 * @param mat_3x2x_view_21 :
	 * - matrix view on \f$M\f$, starting at \f$(2 n_x,n_x)\f$, ending at \f$(3 n_x - 1, 2n_x - 1)\f$ if n > 0
	 * - matrix view on \f$M\f$, starting at \f$(n_x + n_t,n_x)\f$, ending at \f$(2n_x + n_t - 1, n_x + n_t - 1)\f$ else
	 * 
	 * @param vect_2x :
	 * - allocated \f$(2 n_x)\f$-vector if n > 0
	 * - allocated \f$(n_x + n_t)\f$-vector else 
	 * @brief
	 This function computes the square root of the current smoothing state covariance and \f$ [P_{n + 1|N}]^{\frac{1}{2}}\ \; K_{n|N}^T \f$
	 - 1. Building the matrix
	 if n > 0
	\f$
		\begin{pmatrix}
			[Q_2^{x,x}]^{\frac{1}{2}}				&	0						\\
			[P_{n|n}]^{\frac{1}{2}} [F_2^{x,x}]' 	&	[P_{n|n}]^{\frac{1}{2}} \\
			0										&	[P_{n+1|N}]^{\frac{1}{2}} K_{n|N}' 
		\end{pmatrix}
	\f$ 
	else
	\f$
		\begin{pmatrix}
			[Q_2^{x,x}]^{\frac{1}{2}}				&	0						\\
			[Q_{0|0}]^{\frac{1}{2}} [F_2^{x,t}]' 	&	[Q_{0|0}]^{\frac{1}{2}} \\
			0										&	[P_{n+1|N}]^{\frac{1}{2}} K_{n|N}'
		\end{pmatrix}
	\f$ 
	 - 2.QR decomposition of the previous matrix which gives
	 if n > 0
	\f$
		\begin{pmatrix}
			*	&	[P_{n + 1|N}]^{\frac{1}{2}}\ \; K_{n|N}^T		\\
			0 	&	[P_{n|N}]^{\frac{1}{2}} \\
			0	&	0
		\end{pmatrix}
	\f$ 
	else
	\f$
		\begin{pmatrix}
			*	&	[P_{n + 1|N}]^{\frac{1}{2}}\ \; K_{n|N}^T		\\
			0 	&	[Q_{0|N}]^{\frac{1}{2}} \\
			0	&	0
		\end{pmatrix}
	\f$ 
	 * 
	 * 
	**/
	void tkalman_get_sqrt_p_s_and_c_s(	gsl_matrix * sqrt_p_s,
										gsl_matrix * c_s,
										const gsl_matrix * sqrt_p_f,
										const gsl_matrix * sqrt_p_s_,
										const gsl_matrix * f2_xx,
										const gsl_matrix * sqrt_q2_xx,
										const gsl_matrix * s_gain,
										gsl_matrix * mat_3x2x,
										gsl_matrix * mat_3x2x_view_00,
										gsl_matrix * mat_3x2x_view_01,
										gsl_matrix * mat_3x2x_view_10,
										gsl_matrix * mat_3x2x_view_11,
										gsl_matrix * mat_3x2x_view_20,
										gsl_matrix * mat_3x2x_view_21,
										gsl_vector * vect_2x);
	/**@fn void tkalman_do_smoothing(gsl_vector * x_s,
								  gsl_matrix * sqrt_p_s,
								  gsl_matrix * c_s,
								  const gsl_vector * x_f,
								  const gsl_matrix * sqrt_p_f,
								  const gsl_vector * x_p_,
								  const gsl_matrix * sqrt_p_p_,
								  const gsl_vector * x_s_,
								  const gsl_matrix * sqrt_p_s_,
								  const gsl_matrix * f2_xx,
								  const gsl_matrix * sqrt_q2_xx,
								  gsl_matrix * mat_xx,
								  gsl_matrix * mat_3x2x,
								  gsl_matrix * mat_3x2x_view_00,
								  gsl_matrix * mat_3x2x_view_01,
								  gsl_matrix * mat_3x2x_view_10,
								  gsl_matrix * mat_3x2x_view_11,
								  gsl_matrix * mat_3x2x_view_20,
								  gsl_matrix * mat_3x2x_view_21,
								  gsl_permutation * perm_x,
								  gsl_vector * vect_2x)
	 * @param[out] x_s : 
	 \f$ \hat{x}_{n|N} \f$, smoothing state expectation
	 * @param[out] sqrt_p_s :
	 \f$ [P_{n|N}]^{\frac{1}{2}}\f$, square root of the current smoothing state covariance
	 * @param[out] c_s : 
	 \f$ [P_{n + 1|N}]^{\frac{1}{2}}\ \; K_{n|N}^T \f$
	 * @param[in] x_f : 
	 \f$ \hat{x}_{n|n} \f$, filtering state expectation
	 * @param[in] sqrt_p_f : 
	 \f$[P_{n|n}]^{\frac{1}{2}}\f$, square root of the current filtering state covariance matrix
	 * @param[in] x_p_ : 
	 \f$ \hat{x}_{n + 1|n} \f$, following predicting state expectation
	 * @param[in] sqrt_p_p_ : 
	 \f$[P_{n+1|n}]^{\frac{1}{2}}\f$, square root of the following predicting state covariance matrix
	 * @param[in] x_s_ : 
	 \f$ \hat{x}_{n + 1|N} \f$, following smoothing state expectation
	 * @param[in]  sqrt_p_s_ : 
	 \f$ [P_{n + 1|N}]^{\frac{1}{2}}  \f$, square root of the following smoothing state covariance
	 * @param[in] f2_xx :
	 \f$ F_2^{x,x} \f$, \f$ F^{x,x} - Q^{x,y} \: [Q^{y,y}]^{-1} \: F^{y,x} \f$
	 * @param[in] sqrt_q2_xx : 
	 \f$ [Q_2^{x,x}]^{\frac{1}{2}}\f$, square root of the reduced noise covariance.
	 * @param mat_3x2x : 
	 - \f$M\f$, allocated \f$(3 n_x, 2 n_x)\f$-matrix
	 * @param mat_3x2x_view_00 : 
	 * matrix view on \f$M\f$, starting at \f$(0,0)\f$, ending at \f$(n_x - 1, n_x - 1)\f$
	 * @param mat_3x2x_view_01 : 
	 * - matrix view on \f$M\f$, starting at \f$(0,n_x)\f$, ending at \f$(n_x - 1, 2n_x - 1)\f$
	 * @param mat_3x2x_view_10 : 
	 * - matrix view on \f$M\f$, starting at \f$(n_x,0)\f$, ending at \f$(2 n_x - 1, n_x - 1)\f$ if
	 * @param mat_3x2x_view_11 : 
	 * - matrix view on \f$M\f$, starting at \f$(n_x,n_x)\f$, ending at \f$(2 n_x - 1, 2 n_x - 1)\f$
	 * @param mat_3x2x_view_20 :
	 * - matrix view on \f$M\f$, starting at \f$(2 n_x,0)\f$, ending at \f$(3 n_x - 1, n_x - 1)\f$
	 * @param mat_3x2x_view_21 :
	 * - matrix view on \f$M\f$, starting at \f$(2 n_x,n_x)\f$, ending at \f$(3 n_x - 1, 2n_x - 1)\f$
	 * @param perm_x : allocated \f$n_x\f$-permutation
	 * @param vect_2x :allocated \f$(2 n_x)\f$-vector.
	 * @brief
	 This function does the smoothing step of Pairwise Kalman filter.
	 - 1.This function computes Pairwise Kalman smoothing gain:
	  \f$K_{n|N} = P_{n|n} [F_2^{x,x}]' [P_{n+1|n}]^{-1}\f$
	 - 2. Computing the smoothing state expectation:
	  \f$  \hat{x}_{n|N} = \hat{x}_{n|n} + K_{n|N} \: (\hat{x}_{n + 1|N} - \hat{x}_{n+1|n} \f$ 
	 - 3. Computing the square root of the current smoothing state covariance and \f$ [P_{n + 1|N}]^{\frac{1}{2}}\ \; K_{n|N}^T \f$
	 - A. Building the matrix
	\f$
		\begin{pmatrix}
			[Q_2^{x,x}]^{\frac{1}{2}}				&	0						\\
			[P_{n|n}]^{\frac{1}{2}} [F_2^{x,x}]' 	&	[P_{n|n}]^{\frac{1}{2}} \\
			0										&	[P_{n+1|N}]^{\frac{1}{2}} K_{n|N}' 
		\end{pmatrix}
	\f$
	 - B.QR decomposition of the previous matrix which gives
	\f$
		\begin{pmatrix}
			*	&	[P_{n + 1|N}]^{\frac{1}{2}}\ \; K_{n|N}^T		\\
			0 	&	[P_{n|N}]^{\frac{1}{2}} \\
			0	&	0
		\end{pmatrix}
	\f$ 
	 *  
	**/
	void tkalman_do_smoothing(gsl_vector * x_s,
							  gsl_matrix * sqrt_p_s,
							  gsl_matrix * c_s,
							  const gsl_vector * x_f,
							  const gsl_matrix * sqrt_p_f,
							  const gsl_vector * x_p_,
							  const gsl_matrix * sqrt_p_p_,
							  const gsl_vector * x_s_,
							  const gsl_matrix * sqrt_p_s_,
							  const gsl_matrix * f2_xx,
							  const gsl_matrix * sqrt_q2_xx,
							  gsl_matrix * mat_xx,
							  gsl_matrix * mat_3x2x,
							  gsl_matrix * mat_3x2x_view_00,
							  gsl_matrix * mat_3x2x_view_01,
							  gsl_matrix * mat_3x2x_view_10,
							  gsl_matrix * mat_3x2x_view_11,
							  gsl_matrix * mat_3x2x_view_20,
							  gsl_matrix * mat_3x2x_view_21,
							  gsl_permutation * perm_x,
							  gsl_vector * vect_2x);
							  
	/**@fn void tkalman_do_smoothing_0(gsl_vector * t_s,
										gsl_matrix * sqrt_q_s,
										gsl_matrix * c_s,
										const gsl_vector * t_f,
										const gsl_matrix * sqrt_q_f,
										const gsl_vector * x_p_,
										const gsl_matrix * sqrt_p_p_,
										const gsl_vector * x_s_,
										const gsl_matrix * sqrt_p_s_,
										const gsl_matrix * f2_xt,
										const gsl_matrix * sqrt_q2_xx,
										gsl_matrix * mat_tx,
										gsl_matrix * mat_2xpt_xpt,
										gsl_matrix * mat_2xpt_xpt_view_00,
										gsl_matrix * mat_2xpt_xpt_view_01,
										gsl_matrix * mat_2xpt_xpt_view_10,
										gsl_matrix * mat_2xpt_xpt_view_11,
										gsl_matrix * mat_2xpt_xpt_view_20,
										gsl_matrix * mat_2xpt_xpt_view_21,
										gsl_permutation * perm_x,
										gsl_vector * vect_xpt);
	 * @param[out] t_s : 
	 \f$ \hat{t}_{0|N} \f$, smoothing state expectation
	 * @param[out] sqrt_q_s :
	 \f$ [Q_{0|N}]^{\frac{1}{2}}\f$
	 * @param[out] c_s : 
	 \f$ [P_{1|N}]^{\frac{1}{2}}\ \; K_{0|N}^T \f$
	 * @param[in] t_f : 
	 \f$ \hat{t}_{0|0} \f$, filtering state expectation
	 * @param[in] sqrt_q_f : 
	 \f$[Q_{0|0}]^{\frac{1}{2}}\f$, square root of the current filtering state covariance matrix 
	 * @param[in] x_p_ : 
	 \f$ \hat{x}_{1|0} \f$, following predicting state expectation
	 * @param[in] sqrt_p_p_ : 
	 \f$[P_{1|0}]^{\frac{1}{2}}\f$, square root of the following predicting state covariance matrix
	 * @param[in] x_s_ : 
	 \f$ \hat{x}_{1|N} \f$, following smoothing state expectation
	 * @param[in]  sqrt_p_s_ : 
	 \f$ [P_{1|N}]^{\frac{1}{2}}  \f$, square root of the following smoothing state covariance
	 * @param[in] f2_xt :
	 \f$ F_2^{x,t} \f$, \f$ F^{x,t} - Q^{x,y} \: [Q^{y,y}]^{-1} \: F^{y,t} \f$ 
	 * @param[in] sqrt_q2_xx : 
	 \f$ [Q_2^{x,x}]^{\frac{1}{2}}\f$, square root of the reduced noise covariance.
	 * @param mat_tx : 
	 allocated \f$(n_t, n_x)\f$-matrix
	 * @param mat_2xpt_xpt : 
	 - \f$M\f$, allocated \f$(2 n_x + n_t, n_x + n_t)\f$-matrix
	 * @param mat_2xpt_xpt_view_00 : 
	 * matrix view on \f$M\f$, starting at \f$(0,0)\f$, ending at \f$(n_x - 1, n_x - 1)\f$
	 * @param mat_2xpt_xpt_view_01 : 
	 * - matrix view on \f$M\f$, starting at \f$(0,n_x)\f$, ending at \f$(n_x - 1, n_x + n_t - 1)\f$
	 * @param mat_2xpt_xpt_view_10 : 
	 * - matrix view on \f$M\f$, starting at \f$(n_x,0)\f$, ending at \f$(n_x + n_t - 1, n_x - 1)\f$
	 * @param mat_2xpt_xpt_view_11 : 
	 * - matrix view on \f$M\f$, starting at \f$(n_x,n_x)\f$, ending at \f$(n_x + n_t - 1, n_x + n_t - 1)\f$
	 * @param mat_2xpt_xpt_view_20 :
	 * - matrix view on \f$M\f$, starting at \f$(n_x + n_t,0)\f$, ending at \f$(2n_x + n_t - 1, n_x - 1)\f$
	 * @param mat_2xpt_xpt_view_21 :
	 * - matrix view on \f$M\f$, starting at \f$(n_x + n_t,n_x)\f$, ending at \f$(2n_x + n_t - 1, n_x + n_t - 1)\f$
	 * @param perm_x : allocated \f$n_x\f$-permutation
	 * @param vect_2x :allocated \f$(2 n_x)\f$-vector.
	 * @brief
	  This function does the smoothing step of Pairwise Kalman filter.
	 - 1.This function computes Pairwise Kalman smoothing gain:
	  \f$K_{0|N} = Q_{0|0} [F_2^{x,t}]' [P_{1|0}]^{-1}\f$
	 - 2. Computing the smoothing state expectation:
	  \f$  \hat{t}_{0|N} = \hat{x}_{0|0} + K_{0|N} \: (\hat{x}_{1|N} - \hat{x}_{1|0} \f$ 
	 - 3. Computing the square root of the current smoothing state covariance and \f$ [P_{1|N}]^{\frac{1}{2}}\ \; K_{0|N}^T \f$
	 - A. Building the matrix
	\f$
		\begin{pmatrix}
			[Q_2^{x,x}]^{\frac{1}{2}}				&	0						\\
			[Q_{0|0}]^{\frac{1}{2}} [F_2^{x,t}]' 	&	[Q_{0|0}]^{\frac{1}{2}} \\
			0										&	[P_{1|N}]^{\frac{1}{2}} K_{0|N}'
		\end{pmatrix}
	\f$
	 - B.QR decomposition of the previous matrix which gives
	\f$
		\begin{pmatrix}
			*	&	[P_{1|N}]^{\frac{1}{2}}\ \; K_{0|N}^T		\\
			0 	&	[P_{0|N}]^{\frac{1}{2}} \\
			0	&	0
		\end{pmatrix}
	\f$ 
	**/
	void tkalman_do_smoothing_0(gsl_vector * t_s,
								gsl_matrix * sqrt_q_s,
								gsl_matrix * c_s,
								const gsl_vector * t_f,
								const gsl_matrix * sqrt_q_f,
								const gsl_vector * x_p_,
								const gsl_matrix * sqrt_p_p_,
								const gsl_vector * x_s_,
								const gsl_matrix * sqrt_p_s_,
								const gsl_matrix * f2_xt,
								const gsl_matrix * sqrt_q2_xx,
								gsl_matrix * mat_tx,
								gsl_matrix * mat_2xpt_xpt,
								gsl_matrix * mat_2xpt_xpt_view_00,
								gsl_matrix * mat_2xpt_xpt_view_01,
								gsl_matrix * mat_2xpt_xpt_view_10,
								gsl_matrix * mat_2xpt_xpt_view_11,
								gsl_matrix * mat_2xpt_xpt_view_20,
								gsl_matrix * mat_2xpt_xpt_view_21,
								gsl_permutation * perm_x,
								gsl_vector * vect_xpt);
	/**@class tkalman_smoothing
	 * @author 
     * Valérian Némesin
	 * @brief
	 Pairwise Kalman smoothing
	 */
	class tkalman_smoothing
	{
		public:
			/**@fn tkalman_smoothing :: tkalman_smoothing(const gsl_matrix * f2_xt,
														    const gsl_matrix * sqrt_q2_xx)
			 * @param[in] f2_xt :
			 \f$ F_2^{x,t} \f$, \f$ F^{x,t} - Q^{x,y} \: [Q^{y,y}]^{-1} \: F^{y,t} \f$ 
			 * @param[in] sqrt_q2_xx : 
			 \f$ [Q_2^{x,x}]^{\frac{1}{2}}\f$, square root of the reduced noise covariance.
			 */
			tkalman_smoothing(const gsl_matrix * f2_xt,
							   const gsl_matrix * sqrt_q2_xx);
		
			/**@fn int tkalman_smoothing :: setup(const gsl_matrix * f2_xt,
												   const gsl_matrix * sqrt_q2_xx)
			 * @param[in] f2_xt :
			 \f$ F_2^{x,t} \f$, \f$ F^{x,t} - Q^{x,y} \: [Q^{y,y}]^{-1} \: F^{y,t} \f$ 
			 * @param[in] sqrt_q2_xx : 
			 \f$ [Q_2^{x,x}]^{\frac{1}{2}}\f$, square root of the reduced noise covariance.
			 * @return
			 * - 0 OK
			 * - 1 problem
			 * @brief
			 * Setup
			**/
			virtual int setup(const gsl_matrix * f2_xt,
							  const gsl_matrix * sqrt_q2_xx);
		
			/**@fn tkalman_smoothing :: ~ tkalman_smoothing()
			 */
			virtual ~tkalman_smoothing();
			
			/**@fn bool tkalman_smoothing :: operator!() const;
			 * @return
			 * - 0 OK
			 * - 1 else
			 */
			virtual bool operator!() const;
			
			
			/**@fn inline void tkalman_smoothing :: do_smoothing (	gsl_vector * x_s,
																	gsl_matrix * sqrt_p_s,
																	gsl_matrix * c_s,
																	const gsl_vector * x_f,
																	const gsl_matrix * sqrt_p_f,
																	const gsl_vector * x_p_,
																	const gsl_matrix * sqrt_p_p_,
																	const gsl_vector * x_s_,
																	const gsl_matrix * sqrt_p_s_)
			 * @param[out] x_s : 
			 \f$ \hat{x}_{n|N} \f$, smoothing state expectation
			 * @param[out] sqrt_p_s :
			 \f$ [P_{n|N}]^{\frac{1}{2}}\f$, square root of the current smoothing state covariance
			 * @param[out] c_s : 
			 \f$ [P_{n + 1|N}]^{\frac{1}{2}}\ \; K_{n|N}^T \f$
			 * @param[in] x_f : 
			 \f$ \hat{x}_{n|n} \f$, filtering state expectation
			 * @param[in] sqrt_p_f : 
			 \f$[P_{n|n}]^{\frac{1}{2}}\f$, square root of the current filtering state covariance matrix
			 * @param[in] x_p_ : 
			 \f$ \hat{x}_{n + 1|n} \f$, following predicting state expectation
			 * @param[in] sqrt_p_p_ : 
			 \f$[P_{n+1|n}]^{\frac{1}{2}}\f$, square root of the following predicting state covariance matrix
			 * @param[in] x_s_ : 
			 \f$ \hat{x}_{n + 1|N} \f$, following smoothing state expectation
			 * @param[in]  sqrt_p_s_ : 
			 \f$ [P_{n + 1|N}]^{\frac{1}{2}}  \f$, square root of the following smoothing state covariance
			 * @brief
			 This function does the smoothing step of Pairwise Kalman filter.
			 - 1.This function computes Pairwise Kalman smoothing gain:
			  \f$K_{n|N} = P_{n|n} [F_2^{x,x}]' [P_{n+1|n}]^{-1}\f$
			 - 2. Computing the smoothing state expectation:
			  \f$  \hat{x}_{n|N} = \hat{x}_{n|n} + K_{n|N} \: (\hat{x}_{n + 1|N} - \hat{x}_{n+1|n} \f$ 
			 - 3. Computing the square root of the current smoothing state covariance and \f$ [P_{n + 1|N}]^{\frac{1}{2}}\ \; K_{n|N}^T \f$
			 - A. Building the matrix
			\f$
				\begin{pmatrix}
					[Q_2^{x,x}]^{\frac{1}{2}}				&	0						\\
					[P_{n|n}]^{\frac{1}{2}} [F_2^{x,x}]' 	&	[P_{n|n}]^{\frac{1}{2}} \\
					0										&	[P_{n+1|N}]^{\frac{1}{2}} K_{n|N}' 
				\end{pmatrix}
			\f$
			 - B.QR decomposition of the previous matrix which gives
			\f$
				\begin{pmatrix}
					*	&	[P_{n + 1|N}]^{\frac{1}{2}}\ \; K_{n|N}^T		\\
					0 	&	[P_{n|N}]^{\frac{1}{2}} \\
					0	&	0
				\end{pmatrix}
			\f$ 
			 */
			inline void do_smoothing (gsl_vector * x_s,
									  gsl_matrix * sqrt_p_s,
									  gsl_matrix * c_s,
									  const gsl_vector * x_f,
									  const gsl_matrix * sqrt_p_f,
									  const gsl_vector * x_p_,
									  const gsl_matrix * sqrt_p_p_,
									  const gsl_vector * x_s_,
									  const gsl_matrix * sqrt_p_s_)
			{
				tkalman_do_smoothing(x_s,
									 sqrt_p_s,
									 c_s,
								     x_f,
									 sqrt_p_f,
								     x_p_,
									 sqrt_p_p_,
									 x_s_,
									 sqrt_p_s_,
								     &f2_xx,
									 _sqrt_q2_xx,
								     &mat_xx,
							         &mat_3x2x,
							         &mat_3x2x_view_00,
							         &mat_3x2x_view_01,
							         &mat_3x2x_view_10,
							         &mat_3x2x_view_11,
							         &mat_3x2x_view_20,
							         &mat_3x2x_view_21,
								     perm_x,
							         &vect_2x);
			}
								 
			/**@fn inline void tkalman_smoothing :: do_smoothing_0 ( gsl_vector * t_s,
																	 gsl_matrix * sqrt_q_s,
																	 gsl_matrix * c_s,
																	 const gsl_vector * t_f,
																	 const gsl_matrix * sqrt_q_f,
																	 const gsl_vector * x_p_,
																	 const gsl_matrix * sqrt_p_p_,
																	 const gsl_vector * x_s_,
																	 const gsl_matrix * sqrt_p_s_)
			 * @param[out] t_s : 
			 \f$ \hat{t}_{0|N} \f$, smoothing state expectation
			 * @param[out] sqrt_q_s :
			 \f$ [Q_{0|N}]^{\frac{1}{2}}\f$
			 * @param[out] c_s : 
			 \f$ [P_{1|N}]^{\frac{1}{2}}\ \; K_{0|N}^T \f$
			 * @param[in] t_f : 
			 \f$ \hat{t}_{0|0} \f$, filtering state expectation
			 * @param[in] sqrt_q_f : 
			 \f$[Q_{0|0}]^{\frac{1}{2}}\f$, square root of the current filtering state covariance matrix 
			 * @param[in] x_p_ : 
			 \f$ \hat{x}_{1|0} \f$, following predicting state expectation
			 * @param[in] sqrt_p_p_ : 
			 \f$[P_{1|0}]^{\frac{1}{2}}\f$, square root of the following predicting state covariance matrix
			 * @param[in] x_s_ : 
			 \f$ \hat{x}_{1|N} \f$, following smoothing state expectation
			 * @param[in]  sqrt_p_s_ : 
			 \f$ [P_{1|N}]^{\frac{1}{2}}  \f$, square root of the following smoothing state covariance
			 * @brief
			  This function does the smoothing step of Pairwise Kalman filter.
			 - 1.This function computes Pairwise Kalman smoothing gain:
			  \f$K_{0|N} = Q_{0|0} [F_2^{x,t}]' [P_{1|0}]^{-1}\f$
			 - 2. Computing the smoothing state expectation:
			  \f$  \hat{t}_{0|N} = \hat{x}_{0|0} + K_{0|N} \: (\hat{x}_{1|N} - \hat{x}_{1|0} \f$ 
			 - 3. Computing the square root of the current smoothing state covariance and \f$ [P_{1|N}]^{\frac{1}{2}}\ \; K_{0|N}^T \f$
			 - A. Building the matrix
			\f$
				\begin{pmatrix}
					[Q_2^{x,x}]^{\frac{1}{2}}				&	0						\\
					[Q_{0|0}]^{\frac{1}{2}} [F_2^{x,t}]' 	&	[Q_{0|0}]^{\frac{1}{2}} \\
					0										&	[P_{1|N}]^{\frac{1}{2}} K_{0|N}'
				\end{pmatrix}
			\f$
			 - B.QR decomposition of the previous matrix which gives
			\f$
				\begin{pmatrix}
					*	&	[P_{1|N}]^{\frac{1}{2}}\ \; K_{0|N}^T		\\
					0 	&	[P_{0|N}]^{\frac{1}{2}} \\
					0	&	0
				\end{pmatrix}
			\f$ 
			**/
			inline void do_smoothing_0 ( gsl_vector * t_s,
										 gsl_matrix * sqrt_q_s,
									 	 gsl_matrix * c_s,
										 const gsl_vector * t_f,
										 const gsl_matrix * sqrt_q_f,
										 const gsl_vector * x_p_,
									 	 const gsl_matrix * sqrt_p_p_,
								 		 const gsl_vector * x_s_,
										 const gsl_matrix * sqrt_p_s_)
			{
				tkalman_do_smoothing_0 (t_s,
										sqrt_q_s,
										c_s,
										t_f,
										sqrt_q_f,
										x_p_,
										sqrt_p_p_,
										x_s_,
										sqrt_p_s_,
										_f2_x,
										_sqrt_q2_xx,
										mat_tx,
										mat_2xpt_xpt,
										&mat_2xpt_xpt_view_00,
										&mat_2xpt_xpt_view_01,
										&mat_2xpt_xpt_view_10,
										&mat_2xpt_xpt_view_11,
										&mat_2xpt_xpt_view_20,
										&mat_2xpt_xpt_view_21,
										perm_x,
										vect_xpt);
			}
			
			
		protected:
			/**@fn int tkalman_smoothing :: set_params(	const gsl_matrix * f2_x,
														const gsl_matrix * sqrt_q2_xx);
			 * @param[in] f2_xt :
			 \f$ F_2^{x,t} \f$, \f$ F^{x,t} - Q^{x,y} \: [Q^{y,y}]^{-1} \: F^{y,t} \f$ 
			 * @param[in] sqrt_q2_xx : 
			 \f$ [Q_2^{x,x}]^{\frac{1}{2}}\f$, square root of the reduced noise covariance.
			 */
			int set_params(const gsl_matrix * f2_x,
						   const gsl_matrix * sqrt_q2_xx);
		
		
			/**@fn void tkalman_smoothing :: free();
			 * @brief
			 * This function frees memory.
			 */
			void free();
		
			/**@fn int tkalman_smoothing :: alloc();
			 * @return
			 * 0 OK
			 * @brief
			 * This function allocates object memory.
			 */
			int alloc();
		
			/**@fn void tkalman_smoothing :: create_views();
			 * @brief
			 * This function creates matrix views.
			 */
			void create_views();
		
		
			/**@fn void tkalman_smoothing :: initialize();
			 * @brief
			 * This function sets object attributes to 0.
			 */
			void initialize();
		//Données propres
			unsigned int _size_x;
			unsigned int _size_y;
			unsigned int _size_t;
			gsl_matrix mat_3x2x;
			gsl_matrix mat_3x2x_view_00;
			gsl_matrix mat_3x2x_view_01;
			gsl_matrix mat_3x2x_view_10;
			gsl_matrix mat_3x2x_view_11;
			gsl_matrix mat_3x2x_view_20;
			gsl_matrix mat_3x2x_view_21;
			gsl_vector vect_2x;
			gsl_matrix mat_xx;
			gsl_matrix mat_2xpt_xpt_view_00;
			gsl_matrix mat_2xpt_xpt_view_01;
			gsl_matrix mat_2xpt_xpt_view_10;
			gsl_matrix mat_2xpt_xpt_view_11;
			gsl_matrix mat_2xpt_xpt_view_20;
			gsl_matrix mat_2xpt_xpt_view_21;
			
			
			gsl_permutation * perm_x;
			gsl_vector * vect_xpt;
			gsl_matrix * mat_tx;
			gsl_matrix * mat_2xpt_xpt;
		//Paramètres
			const gsl_matrix * _f2_x;
			gsl_matrix f2_xx;
			gsl_matrix f2_xy;
			const gsl_matrix * _sqrt_q2_xx;
	};




#endif
