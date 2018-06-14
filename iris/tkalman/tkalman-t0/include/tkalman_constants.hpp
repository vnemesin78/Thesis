/**@file tkalman_constants.hpp
 * @author Valérian Némesin
 * @brief
 Pairwise Kalman constants
 */
#ifndef _TKALMAN_CONSTANTS2_HPP_
	#define _TKALMAN_CONSTANTS2_HPP_
	#include <gsl/gsl_matrix.h>
	#include <gsl/gsl_linalg.h>
	#include <gsl/gsl_blas.h>
	#include "gsl_triangle_matrix.hpp"
	
	/**@fn void tkalman_get_f2_x_(	gsl_matrix * f2_x_,
									const gsl_matrix * f_x_,
									const gsl_matrix * f_y_,
									const gsl_matrix * q2_xy)
	 * @param[out] f2_x_ : \f$ F_2^{x,z} \f$, \f$ F^{x,z} - Q^{x,y} [Q^{y,y}]^{-1} F^{y,z} \f$, with z = x or y.
	 * @param[in] f_x_ : \f$ F^{x,z} \f$, submatrix of transition matrix.
	 * @param[in] f_y_ : \f$ F^{y,z} \f$, submatrix of transition matrix.
	 * @param[in] q2_xy : \f$ Q_2^{x,y}\f$, \f$ Q^{x,y} [Q^{y,y}]^{-\frac{1}{2}} \f$
	 * @brief
	 * This function computes \f$ F_2^{x,z} \f$ by the formula:
	 * \f$ Q_2^{x,y}\f$, \f$ Q^{x,y} [Q^{y,y}]^{-\frac{1}{2}} \f$
	 */
	void tkalman_get_f2_x_(gsl_matrix * f2_x_,
						   const gsl_matrix * f_x_,
						   const gsl_matrix * f_y_,
						   const gsl_matrix * q2_xy);

	/**@fn void tkalman_get_sqrt_q2_xx_sqrt_q_yy_and_q2_xy_from_sqrt_q(	gsl_matrix * sqrt_q2_xx,
																		gsl_matrix * q2_xy,
																		gsl_matrix * sqrt_q_yy,
																		const gsl_matrix * sqrt_q_view_xx,
																		const gsl_matrix * sqrt_q_view_xy,
																		const gsl_matrix * sqrt_q_view_yy,
																		gsl_matrix * mat_tt,
																		gsl_matrix * mat_tt_yy,
																		gsl_matrix * mat_tt_yx,
																		gsl_matrix * mat_tt_xy,
																		gsl_matrix * mat_tt_xx,
																		gsl_vector * vect_t,
																		gsl_permutation * perm_y)
	 * @param[out] sqrt_q2_xx : 
	 \f$ [Q_2^{x,x}]^{\frac{1}{2}}\f$, square root of the reduced noise covariance.
	 * @param[out] q2_xy : 
	 \f$ Q_2^{x,y} \f$, \f$ Q^{x,y} \: [Q^{y,y}]^{-\frac{1}{2}} \f$
	 * @param[out] sqrt_q_yy : \f$ [Q^{y,y}]^{\frac{1}{2}}\f$, square root of measurement covariance matrix.
	 * @param[in] sqrt_q_view_xx : 
	 \f$ [Q^{\frac{1}{2}}]^{x,x} \f$, matrix view on \f$[Q]^{\frac{1}{2}}\f$,  starting at \f$ (0,0) \f$, ending at \f$ (n_x -1, n_x-1) \f$
	 * @param[in] sqrt_q_view_xy : 
	 \f$ [Q^{\frac{1}{2}}]^{x,y} \f$, matrix view on \f$[Q]^{\frac{1}{2}}\f$,  starting at \f$ (0,n_x) \f$, ending at \f$ (n_t -1, n_x-1) \f$
	 * @param[in] sqrt_q_view_yy :  
	 \f$ [Q^{\frac{1}{2}}]^{y,y}\f$, matrix view on \f$[Q]^{\frac{1}{2}}\f$,  starting at \f$ (n_x,n_x) \f$, ending at \f$ (n_t -1, n_t -1) \f$
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
	 * @param vect_t : 
	 * allocated \f$ n_t-\f$vector.
	 * @param perm_y : 
	 * allocated \f$ n_y-\f$ permutation.
	 * 
	 * @brief
	 This function computes:
	 - \f$ [Q_2^{x,x}]^{\frac{1}{2}}\f$
	 - \f$ Q_2^{x,y} \f$, \f$ Q^{x,y} \: [Q^{y,y}]^{-\frac{1}{2}} \f$
	 - \f$ [Q^{y,y}]^{\frac{1}{2}}\f$
	 - A. First, it builds the matrix :
	  \f$ 
	  M=
		 \begin{pmatrix}
			[Q^{\frac{1}{2}}]^{y,y}	&	0 \\
			[Q^{\frac{1}{2}}]^{x,y} & 	[Q^{\frac{1}{2}}]^{x,x}
		\end{pmatrix}
	 \f$
	 - B. Its effectuates the QR decompostion of the previous matrix which gives:
	 \f$ 
	 \begin{pmatrix}
		[Q^{y,y}]^{\frac{1}{2}}	&	[Q_2^{x,y}]'			\\
	 	0						& 	[Q_2^{x,x}]^{\frac{1}{2}}	 
	 \end{pmatrix}
	\f$ 
	**/
	void tkalman_get_sqrt_q2_xx_sqrt_q_yy_and_q2_xy_from_sqrt_q(gsl_matrix * sqrt_q2_xx,
																gsl_matrix * q2_xy,
																gsl_matrix * sqrt_q_yy,
																const gsl_matrix * sqrt_q_view_xx,
																const gsl_matrix * sqrt_q_view_xy,
																const gsl_matrix * sqrt_q_view_yy,
																gsl_matrix * mat_tt,
																gsl_matrix * mat_tt_yy,
																gsl_matrix * mat_tt_yx,
																gsl_matrix * mat_tt_xy,
																gsl_matrix * mat_tt_xx,
																gsl_vector * vect_t,
																gsl_permutation * perm_y);
	/**@fn void tkalman_get_constants(gsl_matrix * f2_xt,
									  gsl_matrix * sqrt_q2_xx,
									  gsl_matrix * q2_xy,
									  gsl_matrix * sqrt_q_yy,
									  const gsl_matrix * f_xt,
									  const gsl_matrix * f_yt,
									  const gsl_matrix * sqrt_q_view_xx,
									  const gsl_matrix * sqrt_q_view_xy,
									  const gsl_matrix * sqrt_q_view_yy,
									  gsl_matrix * mat_tt,
									  gsl_matrix * mat_tt_yy,
									  gsl_matrix * mat_tt_yx,
									  gsl_matrix * mat_tt_xy,
									  gsl_matrix * mat_tt_xx,
									  gsl_vector * vect_t,
									  gsl_permutation * perm_y)
	 * @param[out] f2_xt : \f$ F_2^{x,t} \f$, \f$ F^{x,t} - Q^{x,y} \: [Q^{y,y}]^{-1} \: F^{y,t} \f$
	 * @param[out] sqrt_q2_xx : \f$ [Q_2^{x,x}]^{\frac{1}{2}}\f$, square root of the reduced noise covariance.
	 * @param[out] q2_xy : \f$ Q_2^{x,y} \f$, \f$ Q^{x,y} \: [Q^{y,y}]^{-\frac{1}{2}} \f$
	 * @param[out] sqrt_q_yy : \f$ [Q^{y,y}]^{\frac{1}{2}}\f$, square root of measurement covariance matrix.
	 * @param[in] f_xt : \f$ F^{x,t} \f$, submatrix of transition matrix.
	 * @param[in] f_yt : \f$ F^{y,t} \f$, submatrix of transition matrix.
	 * @param[in] sqrt_q_view_xx : \f$ [Q^{\frac{1}{2}}]^{x,x} \f$, matrix view on \f$[Q]^{\frac{1}{2}}\f$,  starting at \f$ (0,0) \f$, ending at \f$ (n_x -1, n_x-1) \f$
	 * @param[in] sqrt_q_view_xy : \f$ [Q^{\frac{1}{2}}]^{x,y} \f$, matrix view on \f$[Q]^{\frac{1}{2}}\f$,  starting at \f$ (0,n_x) \f$, ending at \f$ (n_t -1, n_x-1) \f$
	 * @param[in] sqrt_q_view_yy :  \f$ [Q^{\frac{1}{2}}]^{y,y}\f$, matrix view on \f$[Q]^{\frac{1}{2}}\f$,  starting at \f$ (n_x,n_x) \f$, ending at \f$ (n_t -1, n_t -1) \f$
	 * @param mat_tt : \f$M\f$, allocated \f$ ( n_t, n_t) \f$-matrix
	 * @param mat_tt_yy : \f$M^{y,y}\f$, matrix view on \f$M\f$,  starting at \f$ (0,0) \f$, ending at \f$ (n_y -1, n_y - 1) \f$
	 * @param mat_tt_yx : \f$M^{y,x}\f$, matrix view on \f$M\f$,  starting at \f$ (0, n_y) \f$, ending at \f$ (n_y -1, n_t - 1) \f$
	 * @param mat_tt_xy : \f$M^{x,y}\f$, matrix view on \f$M\f$,  starting at \f$ (n_y, 0) \f$, ending at \f$ (n_t - 1, n_y - 1) \f$
	 * @param mat_tt_xx : \f$M^{x,x}\f$, matrix view on \f$M\f$,  starting at \f$ (n_y, n_y) \f$, ending at \f$ (n_t - 1, n_t - 1) \f$
	 * @param vect_t : allocated \f$ n_t-\f$vector.
	 * @param perm_y : allocated \f$ n_y-\f$permutation.
	 * @brief
	 * This function computes Pairwise Kalman constants.
	 * - \f$ F^{x,t} - Q^{t,y} \: [Q^{y,y}]^{-1} \: F^{y,t} \f$
	 * - \f$ [Q_2^{x,x}]^{\frac{1}{2}}\f$
	 * - \f$ Q_2^{x,y} \f$, \f$ Q^{x,y} \: [Q^{y,y}]^{-\frac{1}{2}} \f$
	 * - \f$ [Q^{y,y}]^{\frac{1}{2}}\f$
	 * 
	**/
	void tkalman_get_constants(gsl_matrix * f2_xt,
							   gsl_matrix * sqrt_q2_xx,
							   gsl_matrix * q2_xy,
							   gsl_matrix * sqrt_q_yy,
							   const gsl_matrix * f_xt,
							   const gsl_matrix * f_yt,
							   const gsl_matrix * sqrt_q_view_xx,
							   const gsl_matrix * sqrt_q_view_xy,
							   const gsl_matrix * sqrt_q_view_yy,
							   gsl_matrix * mat_tt,
							   gsl_matrix * mat_tt_yy,
							   gsl_matrix * mat_tt_yx,
							   gsl_matrix * mat_tt_xy,
							   gsl_matrix * mat_tt_xx,
							   gsl_vector * vect_t,
							   gsl_permutation * perm_y);
	
	
	/**@class
	 * @author 
	 * Valérian Némesin
	 * @brief
	 * Pairwise Kalman constants
	 */
	class tkalman_constants
	{
		public:
			/**@fn tkalman_constants :: tkalman_constants(const gsl_matrix * f_xt,
														  const gsl_matrix * f_yt,
														  const gsl_matrix * sqrt_q)
			 * @param[in] f_xt : \f$ F^{x,t} \f$, submatrix of transition matrix.
			 * @param[in] f_yt : \f$ F^{y,t} \f$, submatrix of transition matrix.
			 * @param[in] sqrt_q : \f$ [Q]^{\frac{1}{2}} \f$, square root of covariance matrix
			 * @brief
			 * Constructor
			 * 
			 */
			tkalman_constants(const gsl_matrix * f_xt,
							  const gsl_matrix * f_yt,
							  const gsl_matrix * sqrt_q);
		
			/**@fn int tkalman_constants :: setup( const gsl_matrix * f_xt,
												   const gsl_matrix * f_yt,
												   const gsl_matrix * sqrt_q);
			 * @param[in] f_xt : \f$ F^{x,t} \f$, submatrix of transition matrix.
			 * @param[in] f_yt : \f$ F^{y,t} \f$, submatrix of transition matrix.
			 * @param[in] sqrt_q : \f$ [Q]^{\frac{1}{2}} \f$, square root of covariance matrix
			 * @return
			 * - 0 if success
			 * - else
			 * @brief
			 * Setup
			**/
			virtual int setup(const gsl_matrix * f_xt,
							  const gsl_matrix * f_yt,
							  const gsl_matrix * sqrt_q);
		
			/**@fn tkalman_constants :: ~ tkalman_constants()
			 */
			virtual ~tkalman_constants();
			
			
			
			/**@fn bool tkalman_constants :: operator!() const;
			 * @return
			 * - 0 if OK
			 * - 1 else
			 */
			virtual bool operator!() const;
			
			
			/**@fn inline void tkalman_constants :: get_constants (gsl_matrix * f2_xt,
																   gsl_matrix * sqrt_q2_xx,
																   gsl_matrix * q2_xy,
																   gsl_matrix * sqrt_q_yy)
			 * @param[out] f2_xt : \f$ F_2^{x,t} \f$, \f$ F^{x,t} - Q^{x,y} \: [Q^{y,y}]^{-1} \: F^{y,t} \f$
			 * @param[out] sqrt_q2_xx : \f$ [Q_2^{x,x}]^{\frac{1}{2}}\f$, square root of the reduced noise covariance.
			 * @param[out] q2_xy : \f$ Q_2^{x,y} \f$, \f$ Q^{x,y} \: [Q^{y,y}]^{-1}} \f$
			 * @param[out] sqrt_q_yy : \f$ [Q^{y,y}]^{\frac{1}{2}}\f$, square root of measurement covariance matrix.
			 * @brief
			 * This function computes Pairwise Kalman constants.
			 * - \f$ F^{x,t} - Q^{t,y} \: [Q^{y,y}]^{-1} \: F^{y,t} \f$
			 * - \f$ [Q_2^{x,x}]^{\frac{1}{2}}\f$
			 * - \f$ Q_2^{x,y} \f$, \f$ Q^{x,y} \: [Q^{y,y}]^{-\frac{1}{2}} \f$
			 * - \f$ [Q^{y,y}]^{\frac{1}{2}}\f$
			 */
			inline void get_constants (gsl_matrix * f2_xt,
									   gsl_matrix * sqrt_q2_xx,
									   gsl_matrix * q2_xy,
									   gsl_matrix * sqrt_q_yy)
			{
				tkalman_get_constants(f2_xt,
									  sqrt_q2_xx,
									  q2_xy,
									  sqrt_q_yy,
									  _f_x_,
									  _f_y_,
									  &sqrt_q_view_xx,
									  &sqrt_q_view_xy,
						              &sqrt_q_view_yy,
									  mat_tt,
									  &mat_tt_yy,
									  &mat_tt_yx,
									  &mat_tt_xy,
									  &mat_tt_xx,
									  vect_t,
									  perm_y);
			}
			
		protected:
			/**@fn int tkalman_constants :: set_params(const gsl_matrix * f_xt,
													   const gsl_matrix * f_yt,
													   const gsl_matrix * sqrt_q);
			 * @param[in] f_xt : \f$ F^{x,t} \f$, submatrix of transition matrix.
			 * @param[in] f_yt : \f$ F^{y,t} \f$, submatrix of transition matrix.
			 * @param[in] sqrt_q : \f$ [Q]^{\frac{1}{2}} \f$, square root of covariance matrix
			 * @brief
			 * This function changes object parameters.
			 */
			int set_params(const gsl_matrix * f_xt,
						   const gsl_matrix * f_yt,
						   const gsl_matrix * sqrt_q);;
		
		
			/**@fn void tkalman_constants :: free();
			 * @brief
			 * This function frees memory.
			 */
			void free();
		
			/**@fn int tkalman_constants :: alloc();
			 * @return
			 * - 0 OK
			 * - 1 memory error
			 * @brief
			 * This function allocates attributes.
			 */
			 int alloc();
		
			/**@fn void tkalman_constants :: create_views();
			 * @brief
			 * This function creates matrix views.
			 */
			void create_views();
		
		
			/**@fn void tkalman_constants :: initialize();
			 * @brief
			 * This function sets object attributes to 0.
			 */
			void initialize();
		//Données propres
			unsigned int _size_x;
			unsigned int _size_y;
			unsigned int _size_t;

			gsl_matrix * mat_tt;
			gsl_matrix mat_tt_yy;
			gsl_matrix mat_tt_yx;
			gsl_matrix mat_tt_xy;
			gsl_matrix mat_tt_xx;
			gsl_vector * vect_t;
			gsl_permutation * perm_y;
			
		//Paramètres
			const gsl_matrix * _f_x_;
			const gsl_matrix * _f_y_;
			gsl_matrix sqrt_q_view_xx;
			gsl_matrix sqrt_q_view_xy;
			gsl_matrix sqrt_q_view_yy;
	};

#endif

