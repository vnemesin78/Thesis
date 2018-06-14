/**@file tkalman_filtering.hpp
 * @author Valérian Némesin
 * @date 10/11/2011
 * @brief
Pairwise Kalman filtering
**/
#ifndef _TKALMAN_FILTERING2_HPP_
	#define _TKALMAN_FILTERING2_HPP_
	#include <gsl/gsl_matrix.h>
	#include <gsl/gsl_linalg.h>
	#include <gsl/gsl_blas.h>
	#include "gsl_triangle_matrix.hpp"
	
	/**@fn void tkalman_do_filtering(gsl_vector * x_f,
									 gsl_matrix * sqrt_p_f,
									 gsl_vector * innovation,
									 gsl_matrix * sqrt_s,
									 const gsl_vector * x_p,
									 const gsl_matrix * sqrt_p_p,
									 const gsl_vector * y,
									 const gsl_vector * _y,
									 const gsl_matrix * f_yx,
									 const gsl_matrix * f_yy,
									 const gsl_matrix * sqrt_q_yy,
									 gsl_matrix * mat_tt,
									 gsl_matrix * mat_tt_yy,
									 gsl_matrix * mat_tt_yx,
									 gsl_matrix * mat_tt_xy,
									 gsl_matrix * mat_tt_xx,
									 gsl_matrix * mat_xy,
									 gsl_permutation * perm_y,
									 gsl_vector * vect_t)
	 * @param[out] x_f : 
	 \f$ \hat{x}_{n|n} \f$, filtering state expectation
	 
	 * @param[out] sqrt_p_f : 
	 \f$ [P_{n|n}]^{\frac{1}{2}} \f$, square root of filtering state covariance
	 
	 * @param[out] innovation : 
	 \f$ \tilde{y}_{n} \f$, innovation expectation
	 
	 * @param[out] sqrt_s : 
	 \f$ [S_{n}]^{\frac{1}{2}}  \f$, square root of innovation covariance
	 
	 * @param[in] x_p : 
	 \f$ \hat{x}_{n | n - 1} \f$, predicting state expectation
	 * @param[in] sqrt_p_p : 
	 \f$ [P_{n|n-1}]^{\frac{1}{2}}  \f$, square root of predicting state covariance
	 
	 * @param[in] y : 
	 \f$y_n\f$, current observation
	 
	 * @param[in] _y :
	 \f$y_{n-1}\f$, previous observation
	 
	 * @param[in] f_yx : 
	 \f$ F^{y,x} \f$, submatrix of transition matrix.
	 
	 * @param[in] f_yy : 
	 \f$ F^{y,y} \f$, submatrix of transition matrix.
	 
	 * @param[in] sqrt_q_yy : 
	 \f$ [Q^{y,y}]^{\frac{1}{2}}\f$, square root of measurement covariance matrix.
	 
	 * @param mat_tt : 
	 \f$M\f$, allocated \f$ ( n_t, n_t) \f$-matrix
	 
	 * @param mat_tt_yy : 
	 \f$M^{y,y}\f$, matrix view on \f$M\f$,  starting at \f$ (0,0) \f$, ending at \f$ (n_y -1, n_y - 1) \f$
	 
	 * @param mat_tt_yx : 
	 \f$M^{y,x}\f$, matrix view on \f$M\f$,  starting at \f$ (0, n_y) \f$, ending at \f$ (n_y -1, n_t - 1) \f$
	 
	 * @param mat_tt_xy : 
	 \f$M^{x,y}\f$, matrix view on \f$M\f$,  starting at \f$ (n_y, 0) \f$, ending at \f$ (n_t - 1, n_y - 1) \f$
	 
	 * @param mat_tt_xx : 
	 \f$M^{x,x}\f$, matrix view on \f$M\f$,  starting at \f$ (n_y, n_y) \f$, ending at \f$ (n_t - 1, n_t - 1) \f$
	 
	 * @param mat_xy :
	 \f$M\f$, allocated \f$ ( n_x, n_y ) \f$-matrix
	 
	 * @param vect_t : 
	 * allocated \f$ n_t-\f$vector.
	 
	 * @param perm_y : 
	 allocated \f$ n_y-\f$permutation.
	 
	 * @brief
	 This function does the Pairwise Kalman filtering.
	 - 1. Computing innovation expectation
	 \f$ \tilde{y}_{n} = y_n - F^{y,x} \hat{x}_{n | n - 1} - F^{y,y} y_{n-1} \f$
	 - 2. Computing the matrix
	 \f$ 
	 M = 
	 \begin{pmatrix}
		[Q^{y,y}]^{\frac{1}{2}} 				& 0	\\
		[P_{n|n-1}]^{\frac{1}{2}} [F^{y,x}]'	& [P_{n|n-1}]^{\frac{1}{2}}
	 \end{pmatrix}
	 \f$
	 - 3. QR decompostion of the previous matrix, which gives:
	 \f$
	 \begin{pmatrix}
		[S_{n}]^{\frac{1}{2}}		& [S_{n}]^{\frac{1}{2}} K_{n,n}^T	\\
		0							& [P_{n|n}]^{\frac{1}{2}}
	 \end{pmatrix}
	 \f$
	 - 4. Calculation of filtering state expectation:
	 \f$
		\hat{x}_{n|n} = \hat{x}_{n|n - 1} + K_{n,n} \tilde{y}_{n}
	 \f$
	 */
	void tkalman_do_filtering(gsl_vector * x_f,
								  gsl_matrix * sqrt_p_f,
								  gsl_vector * innovation,
								  gsl_matrix * sqrt_s,
								  const gsl_vector * x_p,
								  const gsl_matrix * sqrt_p_p,
								  const gsl_vector * y,
								  const gsl_vector * _y,
								  const gsl_matrix * f_yx,
								  const gsl_matrix * f_yy,
								  const gsl_matrix * sqrt_q_yy,
								  gsl_matrix * mat_tt,
								  gsl_matrix * mat_tt_yy,
								  gsl_matrix * mat_tt_yx,
								  gsl_matrix * mat_tt_xy,
								  gsl_matrix * mat_tt_xx,
								  gsl_matrix * mat_xy,
								  gsl_permutation * perm_y,
								  gsl_vector * vect_t);
								  
								
	/**@fn void tkalman_do_filtering_0(gsl_vector * t_f_0,
								gsl_matrix * sqrt_q_f_0,
								gsl_vector * innovation,
								gsl_matrix * sqrt_s_0,
								const gsl_vector * t_0,
								const gsl_vector * x_0,
								const gsl_matrix * sqrt_q_0,
								const gsl_vector * y_0,
								const gsl_vector * y_m1,
								const gsl_matrix * f_yt,
								const gsl_matrix * f_yx,
								const gsl_matrix * f_yy,
								const gsl_matrix * sqrt_q_yy,
								gsl_matrix * mat_tpy_tpy,
								gsl_matrix * mat_tpy_tpy_view_00,
								gsl_matrix * mat_tpy_tpy_view_01,
								gsl_matrix * mat_tpy_tpy_view_10,
								gsl_matrix * mat_tpy_tpy_view_11,
								gsl_matrix * mat_ty,
								gsl_permutation * perm_y,
								gsl_vector * vect_tpy);
	 * @param[out] t_f_0 : 
	 \f$ \hat{t}_{0|0} \f$, initial filtering state expectation
	 * @param[out] sqrt_q_f_0, : 
	 \f$ [Q_{0|0}]^{\frac{1}{2}} \f$, square root of initial filtering state covariance
	 
	 * @param[out] innovation : 
	 \f$ \tilde{y}_{0} \f$, initial innovation expectation
	 
	 * @param[out] sqrt_s_0 : 
	 \f$ [S_{0}]^{\frac{1}{2}}  \f$, square root of initial innovation covariance
	 
	 * @param[in] t_0, : 
	 \f$ \hat{t}_{0} \f$, initial (predicting) state expectation

	 * @param[in] x_0 : 
	 \f$ \hat{x}_{0}\f$, initial (predicting) state expectation
	 
	 @param[in] sqrt_q_0 :
	 \f$ [Q_{0}]^{\frac{1}{2}}  \f$, square root of initial (predicting) state covariance
	 
	 * @param[in] y_0 : 
	 \f$y_0\f$, initial observation
	 
	 * @param[in] y_m1 : 
	 \f$y_{-1}\f$, observation "-1" expectation
	 
	 * @param[in] f_yt : 
	 \f$ F^{y,t} \f$, submatrix of transition matrix.
	 
	 * @param[in] f_yx : 
	 \f$ F^{y,x} \f$, submatrix of transition matrix.
	 
	 * @param[in] f_yy : 
	 \f$ F^{y,y} \f$, submatrix of transition matrix.
	 
	 * @param[in] sqrt_q_yy : 
	 \f$ [Q^{y,y}]^{\frac{1}{2}}\f$, square root of measurement covariance matrix.
	 
	 * @param mat_tpy_tpy : 
	 \f$M\f$, allocated \f$ ( n_t + n_y, n_t + n_y ) \f$-matrix
	 
	 * @param mat_tpy_tpy_view_00 : 
	 \f$M^{y,y}\f$, matrix view on \f$M\f$,  starting at \f$ (0,0) \f$, ending at \f$ (n_y -1, n_y - 1) \f$
	 
	 * @param mat_tpy_tpy_view_01 :
	 \f$M^{y,t}\f$, matrix view on \f$M\f$,  starting at \f$ (0,n_y) \f$, ending at \f$ (n_y -1, n_y + n_t - 1) \f$
	 
	 * @param mat_tpy_tpy_view_10 :
	 \f$M^{t,y}\f$, matrix view on \f$M\f$,  starting at \f$ (n_y, 0) \f$, ending at \f$ (n_y + n_t - 1, n_y - 1) \f$
	 
	 * @param mat_tpy_tpy_view_11 :
	 \f$M^{t,t}\f$, matrix view on \f$M\f$,  starting at \f$ (n_y, n_y) \f$, ending at \f$ (n_y + n_t - 1, n_y + n_t - 1) \f$
	 
	 * @param mat_ty : 
	 \f$M\f$, allocated \f$ ( n_t, n_y ) \f$-matrix
	 
	 * @param perm_y : 
	 allocated \f$ n_y-\f$permutation

	 * @param vect_tpy : 
	  allocated \f$ (n_y + n_t)-\f$vector
	 * @brief
	 This function does the Pairwise Kalman filtering.
	 - 1. Computing innovation expectation
	 \f$ \tilde{y}_{0} = y_0 - F^{y,t} \hat{t}_{0} - F^{y,y} \hat{y}_{-1} \f$
	 - 2. Computing the matrix
	 \f$ 
	 M = 
	 \begin{pmatrix}
		[Q^{y,y}]^{\frac{1}{2}} 				& 0	\\
		[Q_{0}]^{\frac{1}{2}} [F^{y,t}]'	& [Q_{0}]^{\frac{1}{2}}
	 \end{pmatrix}
	 \f$
	 - 3. QR decompostion of the previous matrix, which gives:
	 \f$
	 \begin{pmatrix}
		[S_{0}]^{\frac{1}{2}}		& [S_{0}]^{\frac{1}{2}} K_{n,n}^T	\\
		0							& [Q_{0|0}]^{\frac{1}{2}}
	 \end{pmatrix}
	 \f$
	 - 4. Calculation of filtering state expectation:
	 \f$
		\hat{x}_{0|0} = \hat{x}_{0} + K_{0,0} \tilde{y}_{0}
	 \f$
	 */
	void tkalman_do_filtering_0(gsl_vector * t_f_0,
								gsl_matrix * sqrt_q_f_0,
								gsl_vector * innovation,
								gsl_matrix * sqrt_s_0,
								const gsl_vector * t_0,
								const gsl_vector * x_0,
								const gsl_matrix * sqrt_q_0,
								const gsl_vector * y_0,
								const gsl_vector * y_m1,
								const gsl_matrix * f_yt,
								const gsl_matrix * f_yx,
								const gsl_matrix * f_yy,
								const gsl_matrix * sqrt_q_yy,
								gsl_matrix * mat_tpy_tpy,
								gsl_matrix * mat_tpy_tpy_view_00,
								gsl_matrix * mat_tpy_tpy_view_01,
								gsl_matrix * mat_tpy_tpy_view_10,
								gsl_matrix * mat_tpy_tpy_view_11,
								gsl_matrix * mat_ty,
								gsl_permutation * perm_y,
								gsl_vector * vect_tpy);
	
	/**@class tkalman_filtering
	 * @author 
     * Valérian Némesin
	 * @brief
	 * Pairwise Kalman filtering
	 */
	class tkalman_filtering
	{
		public:
			/**@fn tkalman_filtering :: tkalman_filtering(	const gsl_matrix * f_yt,
															const gsl_matrix * sqrt_q_yy);
			 * @param[in] f_yt : 
			 \f$ F^{y,t} \f$, submatrix of transition matrix.
			 * @param[in] sqrt_q_yy : 
			 \f$ [Q^{y,y}]^{\frac{1}{2}}\f$, square root of measurement covariance matrix.
			 */
			tkalman_filtering(	const gsl_matrix * f_yt,
								const gsl_matrix * sqrt_q_yy);
		
			/**@fn int tkalman_filtering :: setup(	const gsl_matrix * f_yt,
													const gsl_matrix * sqrt_q_yy);
			 * @param[in] f_yt : 
			 \f$ F^{y,t} \f$, submatrix of transition matrix.
			 * @param[in] sqrt_q_yy : 
			 \f$ [Q^{y,y}]^{\frac{1}{2}}\f$, square root of measurement covariance matrix.
			 * @return
			 * - 0 OK
			 * - 1 problem
			 * @brief
			 * Setup
			**/
			virtual int setup(	const gsl_matrix * f_yt,
								const gsl_matrix * sqrt_q_yy);
		
			/**@fn tkalman_filtering :: ~ tkalman_filtering()
			 */
			virtual ~tkalman_filtering();
			
			/**@fn bool tkalman_filtering :: operator!() const;
			 * @return
			 * - 0 OK
			 * - 1 memory error
			 * @brief
			 * This method checks object memory.
			 */
			virtual bool operator!() const;
			
			
			/**@fn inline void tkalman_filtering :: do_filtering ( gsl_vector * x_f,
																   gsl_matrix * sqrt_p_f,
																   gsl_vector * innovation,
																   gsl_matrix * sqrt_s,
																   const gsl_vector * x_p,
																   const gsl_matrix * sqrt_p_p,
																   const gsl_vector * y,
																   const gsl_vector * _y)
			 * @param[out] x_f : 
			 \f$ \hat{x}_{n|n} \f$, filtering state expectation
			 * @param[out] sqrt_p_f : 
			 \f$ [P_{n|n}]^{\frac{1}{2}} \f$, square root of filtering state covariance
			 * @param[out] innovation : 
			 \f$ \tilde{y}_{n} \f$, innovation expectation
			 * @param[out] sqrt_s : 
			 \f$ [S_{n}]^{\frac{1}{2}}  \f$, square root of innovation covariance
			 * @param[in] x_p : 
			 \f$ \hat{x}_{n | n - 1} \f$, predicting state expectation 
			 * @param[in] sqrt_p_p : 
			 \f$ [P_{n|n-1}]^{\frac{1}{2}}  \f$, square root of predicting state covariance
			 * @param[in] y : 
			 \f$y_n\f$, current observation
			 * @param[in] _y :
			 \f$y_{n-1}\f$, previous observation
			 * @brief
			 This function does the Pairwise Kalman filtering.
			 - 1. Computing innovation expectation
			 \f$ \tilde{y}_{n} = y_n - F^{y,x} \hat{x}_{n | n - 1} - F^{y,y} y_{n-1} \f$
			 - 2. Computing the matrix
			 \f$ 
			 M = 
			 \begin{pmatrix}
				[Q^{y,y}]^{\frac{1}{2}} 				& 0	\\
				[P_{n|n-1}]^{\frac{1}{2}} [F^{y,x}]'	& [P_{n|n-1}]^{\frac{1}{2}}
			 \end{pmatrix}
			 \f$
			 - 3. QR decompostion of the previous matrix, which gives:
			 \f$
			 \begin{pmatrix}
				[S_{n}]^{\frac{1}{2}}		& [S_{n}]^{\frac{1}{2}} K_{n,n}^T	\\
				0							& [P_{n|n}]^{\frac{1}{2}}
			 \end{pmatrix}
			 \f$
			 - 4. Calculation of filtering state expectation:
			 \f$
				\hat{x}_{n|n} = \hat{x}_{n|n - 1} + K_{n,n} \tilde{y}_{n}
			 \f$
			 */
			inline void do_filtering ( gsl_vector * x_f,
									   gsl_matrix * sqrt_p_f,
									   gsl_vector * innovation,
									   gsl_matrix * sqrt_s,
									   const gsl_vector * x_p,
									   const gsl_matrix * sqrt_p_p,
									   const gsl_vector * y,
									   const gsl_vector * _y)
			{
				tkalman_do_filtering(x_f,
									 sqrt_p_f,
									 innovation,
									 sqrt_s,
									 x_p,
								     sqrt_p_p,
								     y,
								     _y,
									 &f_yx,
									 &f_yy,
									 _sqrt_q_yy,
									 &mat_tt,
                                     &mat_tt_view_00,
                                     &mat_tt_view_01,
                                     &mat_tt_view_10,
                                     &mat_tt_view_11,
                                     &mat_xy,
                                     perm_y,
                                     &vect_t);
			}
								 
			/**@fn inline void tkalman_filtering :: do_filtering_0 ( gsl_vector * t_f_0,
																	 gsl_matrix * sqrt_q_f_0,
																	 gsl_vector * innovation,
																	 gsl_matrix * sqrt_s_0,
																	 const gsl_vector * t_0,
																	 const gsl_vector * x_0,
																	 const gsl_matrix * sqrt_q_0,
																	 const gsl_vector * y_0,
																	 const gsl_vector * y_m1)
			 * @param[out] t_f_0 : 
			 \f$ \hat{t}_{0|0} \f$, initial filtering state expectation
			 * @param[out] sqrt_q_f_0, : 
			 \f$ [Q_{0|0}]^{\frac{1}{2}} \f$, square root of initial filtering state covariance
			 * @param[out] innovation : 
			 \f$ \tilde{y}_{0} \f$, initial innovation expectation
			 * @param[out] sqrt_s_0 : 
			 \f$ [S_{0}]^{\frac{1}{2}}  \f$, square root of initial innovation covariance
			 * @param[in] t_0, : 
			 \f$ \hat{t}_{0} \f$, initial (predicting) state expectation
			 * @param[in] x_0 : 
			 \f$ \hat{x}_{0}  \f$, initial (predicting) state expectation
			 @param[in] sqrt_q_0 :
			 \f$ [Q_{0}]^{\frac{1}{2}}  \f$, square root of initial (predicting) state covariance
			 
			 * @param[in] y_0 : 
			 \f$y_0\f$, initial observation
			 
			 * @param[in] y_m1 : 
			 \f$y_{-1}\f$, observation "-1" expectation
			 * @brief
			 This function does the Pairwise Kalman filtering.
			 - 1. Computing innovation expectation
			 \f$ \tilde{y}_{0} = y_0 - F^{y,t} \hat{t}_{0} - F^{y,y} \hat{y}_{-1} \f$
			 - 2. Computing the matrix
			 \f$ 
			 M = 
			 \begin{pmatrix}
				[Q^{y,y}]^{\frac{1}{2}} 			& 	0	\\
				[Q_{0}]^{\frac{1}{2}} [F^{y,t}]'	& 	[Q_{0}]^{\frac{1}{2}}
			 \end{pmatrix}
			 \f$
			 - 3. QR decompostion of the previous matrix, which gives:
			 \f$
			 \begin{pmatrix}
				[S_{0}]^{\frac{1}{2}}		& [S_{0}]^{\frac{1}{2}} K_{n,n}^T	\\
				0							& [Q_{0|0}]^{\frac{1}{2}}
			 \end{pmatrix}
			 \f$
			 - 4. Calculation of filtering state expectation:
			 \f$
				\hat{x}_{0|0} = \hat{x}_{0} + K_{0,0} \tilde{y}_{0}
			 \f$
			**/
			inline void do_filtering_0 ( gsl_vector * t_f_0,
										 gsl_matrix * sqrt_q_f_0,
										 gsl_vector * innovation,
										 gsl_matrix * sqrt_s_0,
										 const gsl_vector * t_0,
										 const gsl_vector * x_0,
										 const gsl_matrix * sqrt_q_0,
										 const gsl_vector * y_0,
										 const gsl_vector * y_m1)
			{
				tkalman_do_filtering_0 (t_f_0,
										sqrt_q_f_0,
										innovation,
										sqrt_s_0,
										t_0,
										x_0,
										sqrt_q_0,
										y_0,
										y_m1,
										_f_y,
										&f_yx,
										&f_yy,
										_sqrt_q_yy,
										mat_tpy_tpy,
										&mat_tpy_tpy_view_00,
										&mat_tpy_tpy_view_01,
										&mat_tpy_tpy_view_10,
										&mat_tpy_tpy_view_11,
										mat_ty,
										perm_y,
										vect_tpy);
			}
			
			
		protected:
			/**@fn int tkalman_filtering :: set_params(	const gsl_matrix * f_y,
														const gsl_matrix * sqrt_q_yy);
			 * @param[in] f_yt : 
			 \f$ F^{y,t} \f$, submatrix of transition matrix.
			 * @param[in] sqrt_q_yy : 
			 \f$ [Q^{y,y}]^{\frac{1}{2}} \f$, square root of measurement covariance matrix.
			 */
			int set_params(const gsl_matrix * f_y,
						   const gsl_matrix * sqrt_q_yy);
		
		
			/**@fn void tkalman_filtering :: free();
			 * @brief
			 * This function frees memory.
			 */
			void free();
		
			/**@fn int tkalman_filtering :: alloc();
			 * @return
			 * 0 si bon déroulement de l'op.
			 * @brief
			 * This function allocates object memory.
			 */
			int alloc();
		
			/**@fn void tkalman_filtering :: create_views();
			 * @brief
			 * This function creates matrix views.
			 */
			void create_views();
		
		
			/**@fn void tkalman_filtering :: initialize();
			 * @brief
			 * This function sets object attributes to 0.
			 */
			void initialize();
		//Données propres
			unsigned int _size_x;
			unsigned int _size_y;
			unsigned int _size_t;
			
			gsl_matrix * mat_tpy_tpy;
			gsl_matrix * mat_ty;
			gsl_permutation * perm_y;
			gsl_vector * vect_tpy;
			
						
			gsl_matrix mat_tt;
			gsl_matrix mat_tt_view_00;
			gsl_matrix mat_tt_view_01;
			gsl_matrix mat_tt_view_10;
			gsl_matrix mat_tt_view_11;
			gsl_matrix mat_tpy_tpy_view_00;
			gsl_matrix mat_tpy_tpy_view_01;
			gsl_matrix mat_tpy_tpy_view_10;
			gsl_matrix mat_tpy_tpy_view_11;
            gsl_matrix mat_xy;
			
			gsl_vector vect_t;
			
			
			
		//Paramètres
			const gsl_matrix * _f_y;
			gsl_matrix f_yx;
			gsl_matrix f_yy;
			const gsl_matrix * _sqrt_q_yy;
	};
#endif
