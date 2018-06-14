/**@file tkalman_sums.hpp
 * @author Valérian Némesin
 * @brief
Pairwise Kalman sums.
 */

#ifndef _TKALMAN_SUMS2_HPP_
	#define _TKALMAN_SUMS2_HPP_
	#include <gsl/gsl_matrix.h>
	#include <gsl/gsl_linalg.h>
	#include <gsl/gsl_blas.h>
	#include "gsl_triangle_matrix.hpp"
	/**@class tkalman_sums
	 * @author 
     * Valérian Némesin
	 * @brief
	 * Pairwise Kalman sums.
	 * 
	 */
	class tkalman_sums
	{
		public :
			/**@fn tkalman_sums :: tkalman_sums (const gsl_matrix * f2_xt,
												 const gsl_matrix * sqrt_q2_xx);
			 * @param[in] f2_xt :
			 \f$ F_2^{x,t} \f$, \f$ F^{x,t} - Q^{x,y} \: [Q^{y,y}]^{-1} \: F^{y,t} \f$ 
			 * @param[in] sqrt_q2_xx : 
			 \f$ [Q_2^{x,x}]^{\frac{1}{2}}\f$, square root of the reduced noise covariance.
			 * @brief
			 * Constructeur
			 */
			tkalman_sums(const gsl_matrix * f2_xt,
						 const gsl_matrix * sqrt_q2_xx);
						 			 
			/**@fn int tkalman_sums :: setup(const gsl_matrix * f2_xt,
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
						 			 
			/**@fn tkalman_sums :: ~ tkalman_sums()
			 */
			virtual ~tkalman_sums();
			
			/**@fn bool tkalman_sums :: operator!() const;
			 * @return
			 * - 0 OK
			 * - 1 else
			 */
			virtual bool operator!() const;
			
			/**@fn void  tkalman_sums :: compute_sums( const gsl_matrix * const * sqrt_p_f,
													   const gsl_vector * const * x_s,
													   const gsl_matrix * const * sqrt_p_s,
													   const gsl_matrix * const * c_s,
													   const gsl_vector * const * y,
													   unsigned int n,
													   const gsl_vector * x_s_0,
													   const gsl_vector * y_m1,
													   const gsl_matrix * sqrt_q_f_0,
													   const gsl_matrix * c_s_0)
			 * @param[in] sqrt_p_f : 
			 \f$[P_{n|n}]^{\frac{1}{2}}\f$, square roots of the filtering state covariance matrix
			 * @param[in] x_s : 
			 \f$ \hat{x}_{n|N} \f$, smoothing state expectations
			 * @param[in] sqrt_p_s : 
			 \f$[P_{n|N}]^{\frac{1}{2}}\f$, square roots of the smoothing state covariance matrix
			 * @param[in] c_s :  
			 \f$ [P_{n + 1|N}]^{\frac{1}{2}}\ \; K_{n|N}^T \f$
			 * @param[in] y : 
			 \f$ y_{n} \f$, observations
			 * @param[in] n : 
			 number of observations.
			 * @param[in] x_s_0 : 
			 \f$ \hat{x}_{0|N} \f$, initial smoothing state expectation
			 * @param[in] y_m1 : 
			 \f$ \hat{y}_{-1|N} \f$, "initial smoothing state expectation"
			 * @param[in] sqrt_q_f_0 : 
			  \f$[Q_{0|N}]^{\frac{1}{2}}\f$, square roots of the initial smoothing state covariance matrix
			 * @param[in] c_s_0 :  
			 \f$ [P_{1|N}]^{\frac{1}{2}}\ \; K_{0|N}^T \f$
			  
			 * @brief
			 * Computing EM sums.
			 * - A. Initialising 
			 * \f$
			 * 	[W_N]^{\frac{1}{2}} = [Corr(t_{N|N}; t_{N + 1|N}) ]^{\frac{1}{2}} =
			 * \f$
			 * - B. For p = N - 1 to 0
			 *   -- Building the matrix
			 * \f$
			 * M_p = 
			 * \begin{pmatrix}
			 * [W_{p+1}]^{\frac{1}{2}}
			 * [Corr(t_{p|N}; t_{p + 1|N}) ]^{\frac{1}{2}}
			 * \end{pmatrix}
			 * 
			 * \f$
			 * -- QR decompositon which gives:
			 * \f$ 
			 * M_p' = 
			 * \begin{pmatrix}
			 * [W_{p}]^{\frac{1}{2}} \\
			 * 0
			 * \end{pmatrix}
			 * \f$
			 * - C. Result gives:
			 * \f$
			 * M_{sums} = W_{0} = 
			 * [
			 * \begin{pmatrix}
			 * \tilde{C}_{0,0} & \tilde{C}_{0,1} \\
			 * \tilde{C}_{1,0} & \tilde{C}_{1,1}
			 * \end{pmatrix}
			 * ]^{\frac{1}{2}}
			 * \f$
			 */
			void compute_sums(const gsl_matrix * const * sqrt_p_f,
							  const gsl_vector * const * x_s,
							  const gsl_matrix * const * sqrt_p_s,
							  const gsl_matrix * const * c_s,
							  const gsl_vector * const * y,
							  unsigned int n,
							  const gsl_vector * x_s_0,
							  const gsl_vector * y_m1,
							  const gsl_matrix * sqrt_q_f_0,
							  const gsl_matrix * c_s_0);

			/**@fn const gsl_matrix * tkalman_sums :: sqrt_c_00();
			 * @return 
			 \f$
			 M_{sums}^{0,0}
			 \f$
			 */
			inline const gsl_matrix * sqrt_c_00() const
			{
				return &mat_4tp1_2t_view_00;
				
			}
			
			/**@fn const gsl_matrix * tkalman_sums :: sqrt_c_01();
			 * @return 
			 \f$
			 M_{sums}^{0,1}
			 \f$
			 */
			inline const gsl_matrix * sqrt_c_01() const
			{
				return &mat_4tp1_2t_view_01;
				
			}
			
			/**@fn const gsl_matrix * tkalman_sums :: sqrt_c_11();
			 * @return 
			 \f$
			 M_{sums}^{1,1}
			 \f$
			 */
			const gsl_matrix * sqrt_c_11() const
			{
				return &mat_4tp1_2t_view_11;
			}

	protected:
			/**@fn int tkalman_sums :: set_params(const gsl_matrix * f2_x_,
												  const gsl_matrix * sqrt_q2_xx);
			 * @param[in] f2_xt :
			 \f$ F_2^{x,t} \f$, \f$ F^{x,t} - Q^{x,y} \: [Q^{y,y}]^{-1} \: F^{y,t} \f$ 
			 * @param[in] sqrt_q2_xx : 
			 \f$ [Q_2^{x,x}]^{\frac{1}{2}}\f$, square root of the reduced noise covariance.
			 */
			int set_params(const gsl_matrix * f2_x_,
						   const gsl_matrix * sqrt_q2_xx);
		
			/**@fn void tkalman_sums :: free();
			 * @brief
			 * This function frees memory.
			 */
			void free();
		
			/**@fn int tkalman_sums :: alloc();
			 * @return
			 * 0 OK
			 * @brief
			 * This function allocates object memory.
			 */
			int alloc();
		
			/**@fn void tkalman_sums :: create_views();
			 * @brief
			 * This function creates matrix views.
			 */
			void create_views();
		
			/**@fn void tkalman_sums :: initialize();
			 * @brief
			 * This function sets object attributes to 0.
			 */
			void initialize();
			

			

			
			
			
			
		//Données propres
			unsigned int _size_x;
			unsigned int _size_y;
			unsigned int _size_t;
			
			gsl_matrix * mat_2xpt_2xpt;
			
				gsl_matrix mat_2xpt_2xpt_view_00;
				gsl_matrix mat_2xpt_2xpt_view_10;
				gsl_matrix mat_2xpt_2xpt_view_11;
				gsl_matrix mat_2xpt_2xpt_view_12;
				gsl_matrix mat_2xpt_2xpt_view_21;
				gsl_matrix mat_2xpt_2xpt_view_22;
			
				gsl_matrix mat_3x3x;
					gsl_matrix mat_3x3x_view_00;
					gsl_matrix mat_3x3x_view_10;
					gsl_matrix mat_3x3x_view_11;
					gsl_matrix mat_3x3x_view_12;
					gsl_matrix mat_3x3x_view_21;
					gsl_matrix mat_3x3x_view_22;

				gsl_matrix * mat_4tp1_2t;
					gsl_matrix mat_4t2t; // Vue de (0,0) à (4 n_t - 1, 2 n_t - 1)
					gsl_matrix mat_4tp1_2t_view_00; // Vue de (0,0) à (n_t - 1, n_t - 1)
					gsl_matrix mat_4tp1_2t_view_01;	// Vue de (0, n_t) à (n_t - 1, 2 n_t - 1)
					gsl_matrix mat_4tp1_2t_view_11; // Vue de (n_t, n_t) à (2 n_t - 1, 2 n_t - 1)
			
					gsl_matrix mat_2tp1_2t;			// Vue de (2 n_t ,0) à (4 n_t, 2 n_t - 1)
						
						gsl_matrix mat_2tp1_2t_view_00; // Vue de (2 n_t ,0) à (2 n_t + n_x - 1, n_x - 1)
						gsl_matrix mat_2tp1_2t_view_00_bis; // Vue de (2 n_t ,0) à (3 n_t - 1, n_t - 1)
						
						gsl_matrix mat_2tp1_2t_view_02;
						gsl_matrix mat_2tp1_2t_view_02_bis;
						
						gsl_matrix mat_2tp1_2t_view_22;
						
						gsl_vector mat_2tp1_2t_view_30;
						
						gsl_vector mat_2tp1_2t_view_31;
						
						gsl_vector mat_2tp1_2t_view_32;
						
						gsl_vector mat_2tp1_2t_view_33;

				gsl_vector * vect_2t;
				gsl_vector * vect_2xpt;
					gsl_vector vect_3x;
			
		//Paramètres
			const gsl_matrix * _f2_x_;
			gsl_matrix f2_xx;
			gsl_matrix f2_xy;
			const gsl_matrix * _sqrt_q2_xx;
	};
			
/**@fn void tkalman_get_sqrt_cov(const gsl_matrix * sqrt_p_f,
								 const gsl_matrix * sqrt_p_s_,
								 const gsl_matrix * c_s,
								 const gsl_matrix * f2_x_,
								 const gsl_matrix * sqrt_q2_xx,
								 gsl_matrix * mat_3x3x,
								 gsl_matrix * mat_3x3x_view_00,
								 gsl_matrix * mat_3x3x_view_10,
								 gsl_matrix * mat_3x3x_view_11,
								 gsl_matrix * mat_3x3x_view_21,
								 gsl_matrix * mat_3x3x_view_22,
								 gsl_vector * vect_3x);
								 
 * @param[in] sqrt_p_f : 
  \f$[P_{n|n}]^{\frac{1}{2}}\f$, square root of the current filtering state covariance matrix if n > 0
  \f$[Q_{0|0}]^{\frac{1}{2}}\f$, square root of the current filtering state covariance matrix else
 * @param[in]  sqrt_p_s_ : 
  \f$ [P_{n + 1|N}]^{\frac{1}{2}}  \f$, square root of the following smoothing state covariance
 * @param[in] c_s : 
  \f$ [P_{n + 1|N}]^{\frac{1}{2}}\ \; K_{n|N}^T \f$
 * @param[in] f2_x- :
 \f$ F_2^{x,x} \f$, \f$ F^{x,x} - Q^{x,y} \: [Q^{y,y}]^{-1} \: F^{y,x} \f$ if n > 0
 \f$ F_2^{x,t} \f$, \f$ F^{x,t} - Q^{x,y} \: [Q^{y,y}]^{-1} \: F^{y,t} \f$ else
 * @param[in] sqrt_q2_xx : 
 \f$ [Q_2^{x,x}]^{\frac{1}{2}}\f$, square root of the reduced noise covariance.
 * @param mat_3x3x :
 - \f$M\f$, allocated \f$(3 n_x, 3 n_x)\f$-matrix if n > 0
 - \f$M\f$, allocated \f$(2 n_x + n_t, 2 n_x + n_t)\f$-matrix if n = 0
 * @param mat_3x3x_view_00 : 
   matrix view on \f$M\f$, starting at \f$(0,0)\f$, ending at \f$(n_x - 1, n_x - 1)\f$
 * @param mat_3x3x_view_10 : 
	- matrix view on \f$M\f$, starting at \f$(n_x,0)\f$, ending at \f$(2 n_x - 1, n_x - 1)\f$ if n > 0
	- matrix view on \f$M\f$, starting at \f$(n_x,0)\f$, ending at \f$(n_x + n_t - 1, n_x - 1)\f$ else
 * @param mat_3x3x_view_11 :
	 - matrix view on \f$M\f$, starting at \f$(n_x,n_x)\f$, ending at \f$(2 n_x - 1, 2 n_x - 1)\f$ if n > 0
	 - matrix view on \f$M\f$, starting at \f$(n_x,n_x)\f$, ending at \f$(n_x + n_t - 1, n_x + n_t - 1)\f$ else
 * @param mat_3x3x_view_21 :
	- matrix view on \f$M\f$, starting at \f$(2 n_x,n_x)\f$, ending at \f$(3 n_x - 1, 2n_x - 1)\f$ if n > 0
	- matrix view on \f$M\f$, starting at \f$(n_x + n_t,n_x)\f$, ending at \f$(2n_x + n_t - 1, n_x + n_t - 1)\f$ else
 * @param mat_3x3x_view_22 :
	- matrix view on \f$M\f$, starting at \f$(2 n_x,2 n_x)\f$, ending at \f$(3 n_x - 1, 3n_x - 1)\f$ if n > 0
	- matrix view on \f$M\f$, starting at \f$(n_x + n_t,n_x + n_t)\f$, ending at \f$(2n_x + n_t - 1, 2 n_x + n_t - 1)\f$ else
 * @param vect_3x : 
	 * - allocated \f$(3 n_x)\f$-vector if n > 0
	 * - allocated \f$(2 n_x + n_t)\f$-vector else 
 * @brief
 * This function computes the square root of the covariance matric of the random vector \f$ (z_{n|N}; x_{n + 1 |N} ) \f$, where \f$( z = x or t )\f$. 
 * - 1. Building the matrix 
 * \f$
 * M= \begin{pmatrix}
 * \left[Q_2^{x,x}\right]^{\frac{1}{2}}							&	0_{n_x, n_t}	&	0_{n_x, n_x}	\\
 * \left[Z_{n|n}\right]^{\frac{1}{2}} \left[F_2^{x} \right]^T	&	\left[Z_{n|n}\right]^{\frac{1}{2}}	&	0_{n_t, n_x}	\\
 * 0_{n_x, n_x}													&	\left[P_{n + 1 |N} K_{n|N}^T\right]^{\frac{1}{2}}	&	\left[P_{n + 1 |N} K_{n|N}^T\right]^{\frac{1}{2}}
 * \end{pmatrix}
 * \f$
 * - 2. QR decomposition which can be written:
 * \f$ 
   \begin{pmatrix}
 * 	M_{0,0}	&	M_{0,1}	&	M_{0,2} \\
 * 	O	    	&		M_{1,1}	&	M_{1,2} \\ 
 * 	O   	 	&		O	       &	M_{2,2}
 * \end{pmatrix}
 * \f$ 
 * -3. Cov. matrix is 
 * \f$ 
 * \begin{pmatrix}
 * 	M_{1,1}	&	M_{1,2} \\ 
 * 	0    	&	M_{2,2}
 * \end{pmatrix}
 * \f$ 
 * @warning
 * If you want to get the result, you have to build a view on the matrix \f$M\f$,  
 * - starting at \f$(0,0)\f$, ending at \f$(2 n_x - 1, 2 n_x - 1)\f$ if n > 0 
 * - or starting at \f$(0,0)\f$, ending at \f$(n_x + n_t - 1, n_x + n_t - 1)\f$ if n = 0
 *
**/
void tkalman_get_sqrt_cov(const gsl_matrix * sqrt_p_f,
						  const gsl_matrix * sqrt_p_s_,
						  const gsl_matrix * c_s,
						  const gsl_matrix * f2_x_,
						  const gsl_matrix * sqrt_q2_xx,
						  gsl_matrix * mat_3x3x,
						  gsl_matrix * mat_3x3x_view_00,
						  gsl_matrix * mat_3x3x_view_10,
						  gsl_matrix * mat_3x3x_view_11,
						  gsl_matrix * mat_3x3x_view_21,
						  gsl_matrix * mat_3x3x_view_22,
						  gsl_vector * vect_3x);

/**@fn void tkalman_get_sqrt_corr( const gsl_vector * x_s,
									const gsl_vector * x_s_,
									const gsl_vector * _y,
									const gsl_vector * y,
									const gsl_matrix * sqrt_cov_view_00,
									const gsl_matrix * sqrt_cov_view_01,
									const gsl_matrix * sqrt_cov_view_11,
									gsl_matrix * mat_2tp1_2t,
									gsl_matrix * mat_2tp1_2t_view_00,
									gsl_matrix * mat_2tp1_2t_view_02, 
									gsl_matrix * mat_2tp1_2t_view_22,
									gsl_vector * mat_2tp1_2t_view_30,
									gsl_vector * mat_2tp1_2t_view_31,
									gsl_vector * mat_2tp1_2t_view_32,
									gsl_vector * mat_2tp1_2t_view_33,
									gsl_vector * vect_2t );
 * @param[in] x_s : 
 - \f$ \hat{x}_{n|N} \f$, smoothing state expectation if n > 0
 - \f$ \hat{t}_{0|N} \f$, smoothing state expectation if n = 0
 * @param[in] x_s_ : 
 \f$ \hat{x}_{n + 1|N} \f$, following smoothing state expectation
 * @param[in] _y :
 \f$y_{n-1}\f$, previous observation
 * @param[in] y : 
 \f$y_n\f$, current observation
 * @param[in] sqrt_cov_view_00 : 
 * - matrix view on \f$ cov(x_{n|N}; x_{n + 1|N})  \f$, starting at \f$(0, 0)\f$, ending at \f$(n_x - 1, n_x - 1)\f$ if n > 0
 * - matrix view on \f$ cov(t_{0|N}; x_{1|N})  \f$, starting at \f$(0, 0)\f$, ending at \f$(n_t - 1, n_t - 1)\f$ if n = 0
 * @param[in] sqrt_cov_view_01 : 
 * - matrix view on \f$ cov(x_{n|N}; x_{n + 1|N})  \f$, starting at \f$(0, n_x)\f$, ending at \f$(n_x - 1, 2 n_x - 1)\f$ if n > 0
 * - matrix view on \f$ cov(t_{0|N}; x_{1|N})  \f$, starting at \f$(0, n_t)\f$, ending at \f$(n_t - 1, n_t  + n_x - 1)\f$ if n = 0
 * @param[in] sqrt_cov_view_11 : 
 * - matrix view on \f$ cov(x_{n|N}; x_{n + 1|N})  \f$, starting at \f$(n_x, n_x)\f$, ending at \f$(2 n_x - 1, 2 n_x - 1)\f$ if n > 0
 * - matrix view on \f$ cov(t_{0|N}; x_{1|N})  \f$, starting at \f$(n_t, n_t)\f$, ending at \f$(n_t + n_x - 1, n_t + n_x - 1)\f$ if n = 0
 * @param mat_2tp1_2t : 
 * \f$M\f$, allocated \f$(2n_t + 1, 2 n_t)\f$-matrix
 * @param mat_2tp1_2t_view_00 : 
 * matrix view on \f$M\f$ :
 * - starting at \f$(0, 0)\f$, ending at \f$(n_x - 1, n_x - 1)\f$ if n > 0
 * - starting at \f$(0, 0)\f$, ending at \f$(n_t - 1, n_t - 1)\f$ if n = 0
 * @param mat_2tp1_2t_view_02 : 
 * matrix view on \f$M\f$ :
 * - starting at \f$(0, n_t)\f$, ending at \f$(n_x - 1, n_t + n_x - 1)\f$ if n > 0
 * - starting at \f$(0, n_t)\f$, ending at \f$(n_t - 1, n_t + n_x - 1)\f$ if n = 0
 * @param mat_2tp1_2t_view_22 :
 * matrix view on \f$M\f$ :
 * - starting at \f$(n_t, n_t)\f$, ending at \f$(n_t + n_x - 1, n_t + n_x - 1)\f$
 * @param mat_2tp1_2t_view_30 : 
 * matrix view on \f$M\f$ (vect. view):
 * - starting at \f$(2 n_t, 0)\f$, ending at \f$(2 n_t, n_x - 1)\f$
 * @param mat_2tp1_2t_view_31 : 
 * matrix view on \f$M\f$ (vect. view):
 * - starting at \f$(2 n_t, n_x)\f$, ending at \f$(2 n_t, n_t - 1)\f$
 * @param mat_2tp1_2t_view_32 : 
 * matrix view on \f$M\f$ (vect. view):
 * - starting at \f$(2 n_t, n_t)\f$, ending at \f$(2 n_t, n_t + n_x - 1)\f$
 * @param mat_2tp1_2t_view_33 : 
 * matrix view on \f$M\f$ (vect. view):
 * - starting at \f$(2 n_t, n_t + n_x)\f$, ending at \f$(2 n_t, 2 n_t -1 )\f$
 * @param vect_2t : 
 * allocated \f$(2n_t)\f$-vector
 * @brief
 * This function performs the computation of the square root of the correlation matrix of the vector \f$(t_{n|N}; t_{n + 1|N})\f$:
 * - A. Building the matrix :
 * \f$
 * \begin{pmatrix}
 * [[cov]^{\frac{1}{2}}]_{0,0}					&	[[cov]^{\frac{1}{2}}]_{0,1}			&	0	\\
 * 0											&	[[cov]^{\frac{1}{2}}]_{1,1}			&	0	\\
 * 0											&	0									&	0	\\
 * (\hat{x}_{n|N}^T, \hat{y}_{n - 1|N}^T)		&	\hat{x}_{n + 1|N}					&	y_{n}
 * \end{pmatrix}
 * \f$
 * - B. QR decompostion : Result can be written :
 * \f$
 * M = 
 * \begin{pmatrix}
 * M_{0,0}	&	M_{0,1}	&	M_{0,2}	&	M_{0,3}	\\
 * 0		&	M_{1,1}	&	M_{1,2}	&	M_{1,3}	\\
 * 0		&	0		&	M_{2,2}	&	M_{2,3}	\\
 * 0		&	0		&	0		&	M_{3,3}	\\
 * 0		&	0		&	0		&	0		
 * \end{pmatrix}
 * \f$
 * - C. Square root of the correlation matrix is given by:
 * \f$
 * \begin{pmatrix}
 * M_{0,0}	&	M_{0,1}	&	M_{0,2}	&	M_{0,3}	\\
 * 0		&	M_{1,1}	&	M_{1,2}	&	M_{1,3}	\\
 * 0		&	0		&	M_{2,2}	&	M_{2,3}	\\
 * 0		&	0		&	0		&	M_{3,3}
 * \end{pmatrix}
 * \f$
 * 
 * @warning
 * If you want to get the result, you have to build a view on the matrix \f$M\f$,  
 * - starting at \f$(0,0)\f$, ending at \f$(2 n_t - 1, 2 n_t - 1)\f$ 
 */
void tkalman_get_sqrt_corr( const gsl_vector * x_s,
						    const gsl_vector * x_s_,
						    const gsl_vector * _y,
						    const gsl_vector * y,
						    const gsl_matrix * sqrt_cov_view_00,
						    const gsl_matrix * sqrt_cov_view_01,
						    const gsl_matrix * sqrt_cov_view_11,
						    gsl_matrix * mat_2tp1_2t,
						    gsl_matrix * mat_2tp1_2t_view_00, // Modif
						    gsl_matrix * mat_2tp1_2t_view_02, 
						    gsl_matrix * mat_2tp1_2t_view_22, // Modif
						    gsl_vector * mat_2tp1_2t_view_30,
						    gsl_vector * mat_2tp1_2t_view_31,
						    gsl_vector * mat_2tp1_2t_view_32,
						    gsl_vector * mat_2tp1_2t_view_33,
						    gsl_vector * vect_2t );
#endif
