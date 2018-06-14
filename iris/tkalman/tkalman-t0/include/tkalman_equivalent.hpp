/**@file tkalman_equivalent.hpp
 * @author Valérian Némesin
 * @brief
 * Pairwise Kalman equivalent models.
 */
#ifndef _TKALMAN_EQUIVALENT2_HPP_
	#define _TKALMAN_EQUIVALENT2_HPP_
	#include <gsl/gsl_matrix.h>
	#include <gsl/gsl_linalg.h>
	#include <gsl/gsl_blas.h>
	#include "gsl_triangle_matrix.hpp"
	
	/**@fn void tkalman_compute_equivalent_f(gsl_matrix * f_eq,
											 const gsl_matrix * f,
											 const gsl_matrix * m,
											 const gsl_matrix * m_inv,
											 gsl_matrix * mat_tt)
	 * @param[out] f_eq : $F_M$, \f$M\f$-equivalent of \f$F\f$ \f$ (MFM^{-1}) \f$
	 * @param[in] f : \f$F\f$, transition matrix
	 * @param[in] m : \f$M\f$, transformation matrix
	 * @param[in] m_inv : \f$M^{-1}\f$, inverse of transformation matrix
	 * @param mat_tt : allocated \f$ ( n_t, n_t) \f$-matrix
	 * @brief
	 * This functions computes the \f$M\f$-equivalent of \f$F\f$
	 * \f$ F_M = MFM^{-1} \f$
	 * 
	 */
	void tkalman_compute_equivalent_f(gsl_matrix * f_eq,
									  const gsl_matrix * f,
									  const gsl_matrix * m,
									  const gsl_matrix * m_inv,
									  gsl_matrix * mat_tt);
									  
	/**@fn void tkalman_compute_equivalent_sqrt_q(gsl_matrix * sqrt_q_eq,
												  const gsl_matrix * sqrt_q,
												  const gsl_matrix * m,
												  gsl_vector * vect_t)
	 * @param[out] sqrt_q_eq : \f$ [Q_M]^{\frac{1}{2}} \f$, square root of the \f$M\f$-equivalent of covariance matrix \f$ Q \f$
	 * @param[in] sqrt_q : \f$ [Q]^{\frac{1}{2}} \f$, square root of covariance matrix
	 * @param[in] m : \f$M\f$, transformation matrix
	 * @param vect_t : allocated \f$ n_t-\f$vector.
	 * @brief
	 * This function calculates the \f$M\f$-equivalent of covariance matrix \f$ Q \f$.
	 * - Computation of product \f$ [Q]^{\frac{1}{2}} M' \f$
	 * - QR decompostion of previous product
	 * - Resulting matrix gives \f$ [Q_M]^{\frac{1}{2}} \f$.
	 * 
	 */
	void tkalman_compute_equivalent_sqrt_q(gsl_matrix * sqrt_q_eq,
										   const gsl_matrix * sqrt_q,
										   const gsl_matrix * m,
										   gsl_vector * vect_t);
										   
	/**@fn void tkalman_compute_equivalent_t_0(gsl_vector * t_0_eq,
											   const gsl_vector * t_0,
											   const gsl_matrix * m)
	 * @param[out] t_0_eq : \f$ \hat{t}_0^M\f$, \f$M\f$-equivalent of initial state expectation
	 * @param[in] t_0 : \f$ \hat{t}_0 \f$, initial state expectation
	 * @param[in] m : \f$M\f$, transformation matrix
	 * @brief
	 * This function computes the \f$M\f$-equivalent initial state expectation.
	 * - \f$ M \hat{t}_0 \f$
	*/
	void tkalman_compute_equivalent_t_0(gsl_vector * t_0_eq,
										const gsl_vector * t_0,
										const gsl_matrix * m);
										
	/**@fn void tkalman_compute_equivalent_expectation(gsl_vector * x,
													   const gsl_vector * _y,
													   const gsl_matrix * m_xx,
													   const gsl_matrix * m_xy,
													   gsl_vector * vect_x)
	 * @param x : 
	 * - in : state expectation
	 * - out : equivalent state expectation
	 * @param[in] _y : 
	 * - previous observation expectation 
	 * - previous observation value if known
	 * @param[in] m_xx : \f$ M^{x,x} \f$, view on the matrix \f$ M \f$ starting at \f$(0,0)\f$, ending at \f$( n_x - 1, n_x - 1 )\f$
	 * @param[in] m_xy : \f$ M^{x,y} \f$, view on the matrix \f$ M \f$ starting at \f$(0,n_x)\f$, ending at \f$( n_x - 1, n_t - 1 )\f$
	 * @param vect_x : allocated \f$ n_x-\f$vector.
	 * @brief
	 * This function computes the equivalent state expectation:
	 \f$
	 \hat{x}^M_{n|r} = M^{x,x} * \hat{x}_{n|r} + M^{x,y} * \hat{y}_{n - 1|r}
	 \f$
	 */
	void tkalman_compute_equivalent_expectation(gsl_vector * x,
												const gsl_vector * _y,
												const gsl_matrix * m_xx,
												const gsl_matrix * m_xy,
												gsl_vector * vect_x);
												
	/**@fn void tkalman_compute_equivalent_sqrt_cov_tt(gsl_matrix * sqrt_cov_tt,
													   const gsl_matrix * m,
													   gsl_matrix * mat_tt,
													   gsl_vector * vect_t)
	 * @param sqrt_cov_tt :
	 * in : \f$[Q_{n|r}]^{\frac{1}{2}}\f$, square root of covariance matrix of vector \f$t_{n|r}\f$
	 * out : \f$[Q^M_{n|r}]^{\frac{1}{2}}\f$, square root of covariance matrix of \f$ M \f$-equivalent vector \f$t^M_{n|r}\f$
	 * @param[in] m : \f$M\f$, transformation matrix
	 * @param mat_tt : allocated \f$ ( n_t, n_t) \f$-matrix
	 * @param vect_t : allocated \f$ n_t-\f$vector.
	 * @brief
	 * This function computes the square root of covariance matrix of \f$ M \f$-equivalent vector \f$t^M_{n|r}\f$
	 * - Calculation of the product : \f$  [Q_{n|r}]^{\frac{1}{2}} M' \f$
	 * - QR decompostion of the previous product
	 * - Result gives \f$[Q^M_{n|r}]^{\frac{1}{2}}\f$.
	 */
	void tkalman_compute_equivalent_sqrt_cov_tt(gsl_matrix * sqrt_cov_tt,
												const gsl_matrix * m,
												gsl_matrix * mat_tt,
												gsl_vector * vect_t);
																			
	/**@fn void tkalman_compute_equivalent_sqrt_p_with_known_y(gsl_matrix * sqrt_p,
															   const gsl_matrix * m_xx,
															   gsl_matrix * mat_xx,
															   gsl_vector * vect_x)
	 * @param sqrt_p :
	 * in : \f$[P_{n|r}]^{\frac{1}{2}}\f$, square root of covariance matrix of vector \f$x_{n|r}\f$
	 * out : \f$[P^M_{n|r}]^{\frac{1}{2}}\f$, square root of covariance matrix of \f$ M \f$-equivalent vector \f$x^M_{n|r}\f$
	 * @param[in] m_xx : \f$ M^{x,x} \f$, view on the matrix \f$ M \f$ starting at \f$(0,0)\f$, ending at \f$( n_x - 1, n_x - 1 )\f$
	 * @param mat_xx : allocated \f$ ( n_x, n_x) \f$-matrix
	 * @param vect_x : allocated \f$ n_x-\f$vector.
	 * @brief
	 * This function computes the square root of covariance matrix of \f$ M \f$-equivalent vector \f$x^M_{n|r}\f$
	 * - Calculation of the product : \f$  [P_{n|r}]^{\frac{1}{2}} M^{x,x}' \f$
	 * - QR decompostion of the previous product
	 * - Result gives \f$[P^M_{n|r}]^{\frac{1}{2}}\f$.
	 * @warning
	 * Previous observation has to be known. 
	 */
	void tkalman_compute_equivalent_sqrt_p_with_known_y(gsl_matrix * sqrt_p,
														const gsl_matrix * m_xx,
														gsl_matrix * mat_xx,
														gsl_vector * vect_x);
	/**@class
	 * @author 
     * Valérian Némesin
	 * @brief
	 * Equivalent systems
	 * 
	 **/
	class tkalman_equivalents
	{
		public:
			/**@fn tkalman_equivalents :: tkalman_equivalents(	const gsl_matrix * m,
																unsigned int size_x);
			 * @param[in] m : \f$M\f$, transformation matrix
			 * @param[in] size_x : \f$n_x\f$, dimension of vector $\f$x\f$.
			 */
			tkalman_equivalents(const gsl_matrix * m,
							    unsigned int size_x);
		
			/**@fn int tkalman_equivalents :: setup(	const gsl_matrix * m,
																unsigned int size_x);
			 * @param[in] m : \f$M\f$, transformation matrix
			 * @param[in] size_x : \f$n_x\f$, dimension of vector $\f$x\f$.
			 * - 0 OK
			 * - 1 problem
			 * @brief
			 * Setup
			**/
			virtual int setup(const gsl_matrix * m,
							  unsigned int size_x);
		
			/**@fn tkalman_equivalents :: ~ tkalman_equivalents()
			 */
			virtual ~tkalman_equivalents();
			
			/**@fn bool tkalman_equivalents :: operator!() const;
			 * @return
			 * - 0 si OK
			 * - 1 problem
			 * @brief
			 * This method checks object memory.
			 */
			virtual bool operator!() const;
			
			/**@fn inline void tkalman_equivalents :: compute_equivalent_f(gsl_matrix * f)
			 * @param f
			 * - in : \f$ F \f$, transition matrix
			 * - out : \f$ F^{M} \f$, \f$M\f$-equivalent transition matrix
			 * @brief 
			 *  This function computes the \f$M\f$-equivalent transition matrix:
			 * \f$ F_M = MFM^{-1} \f$
			 */
			inline void compute_equivalent_f(gsl_matrix * f)
			{
				gsl_matrix_memcpy(mat_tt_bis, f);
				tkalman_compute_equivalent_f(f,
											 mat_tt_bis,
									         m,
									         m_inv,
									         mat_tt);
			}
			
			
			/**@fn void tkalman_equivalents :: compute_equivalent_sqrt_q(gsl_matrix * sqrt_q)
			 * @param sqrt_q :
			 * - in : \f$ [Q]^{\frac{1}{2}} \f$, square root of covariance matrix
			 * - out : \f$ [Q_M]^{\frac{1}{2}} \f$, square root of the \f$M\f$-equivalent of covariance matrix \f$ Q \f$
			 * @brief
			 * This function calculates the \f$M\f$-equivalent of covariance matrix \f$ Q \f$.
			 * - Computation of product \f$ [Q]^{\frac{1}{2}} M'\f$
			 * - QR decompostion of previous product
			 * - Resulting matrix gives \f$ [Q_M]^{\frac{1}{2}} \f$.
			 * 
			 */
			inline void compute_equivalent_sqrt_q(gsl_matrix * sqrt_q)
			{
				gsl_matrix_memcpy(mat_tt, sqrt_q);
				tkalman_compute_equivalent_sqrt_q(sqrt_q,
												  mat_tt,
												  m,
												  vect_t);
				
			}
			
			
			
			/**@fn void tkalman_equivalents :: compute_equivalent_x(	gsl_vector * x,
																		const gsl_vector * _y)
			 * @param x : 
			 * - in : state expectation
			 * - out : equivalent state expectation
			 * @param[in] _y : 
			 * - previous observation expectation 
			 * - previous observation value if known
			 * @brief
			 * This function computes the equivalent state expectation:
			 * \f$
			 * \hat{x}^M_{n|r} = M^{x,x} * \hat{x}_{n|r} + M^{x,y} * \hat{y}_{n - 1|r}
			 * \f$
			 **/
			inline void compute_equivalent_x(gsl_vector * x,
											 const gsl_vector * _y)
			{
				
				tkalman_compute_equivalent_expectation(x,
													   _y,
													   &m_xx,
													   &m_xy,
													   &vect_x);
			}
			
			/**@fn void tkalman_equivalents :: compute_equivalent_sqrt_cov_tt(gsl_matrix * sqrt_cov_tt)
			 * @param sqrt_cov_tt :
			 * in : \f$[Q_{n|r}]^{\frac{1}{2}}\f$, square root of covariance matrix of vector \f$t_{n|r}\f$
			 * out : \f$[Q^M_{n|r}]^{\frac{1}{2}}\f$, square root of covariance matrix of \f$ M \f$-equivalent vector \f$t^M_{n|r}\f$
			 * @brief
			 * This function computes the square root of covariance matrix of \f$ M \f$-equivalent vector \f$t^M_{n|r}\f$
			 * - Calculation of the product : \f$  [Q_{n|r}]^{\frac{1}{2}} M' \f$
			 * - QR decompostion of the previous product
			 * - Result gives \f$[Q^M_{n|r}]^{\frac{1}{2}}\f$.
			 **/
			inline void compute_equivalent_sqrt_cov_tt(gsl_matrix * sqrt_cov_tt)
			{
				
				tkalman_compute_equivalent_sqrt_cov_tt(sqrt_cov_tt,
													   m,
													   mat_tt,
													   vect_t);	
			}
			
			
			
			/**@fn void tkalman_equivalents :: compute_equivalent_sqrt_p_with_known_y(gsl_matrix * sqrt_p)
			 * @param sqrt_p :
			 * in : \f$[P_{n|r}]^{\frac{1}{2}}\f$, square root of covariance matrix of vector \f$x_{n|r}\f$
			 * out : \f$[P^M_{n|r}]^{\frac{1}{2}}\f$, square root of covariance matrix of \f$ M \f$-equivalent vector \f$x^M_{n|r}\f$
			 * @brief
			 * This function computes the square root of covariance matrix of \f$ M \f$-equivalent vector \f$x^M_{n|r}\f$
			 * - Calculation of the product : \f$  [P_{n|r}]^{\frac{1}{2}} M^{x,x}' \f$
			 * - QR decompostion of the previous product
			 * - Result gives \f$[P^M_{n|r}]^{\frac{1}{2}}\f$.
			 * @warning
			 * Previous observation has to be known. 
			 */
			inline void compute_equivalent_sqrt_p_with_known_y(gsl_matrix * sqrt_p)
			{
				tkalman_compute_equivalent_sqrt_p_with_known_y(sqrt_p,
															   &m_xx,
														       &mat_xx,
														       &vect_x);
			}
			
			
		protected:
			/**@fn int tkalman_equivalents :: set_params(	const gsl_matrix * m,
															unsigned int size_x);
			 * @param[in] m : \f$M\f$, transformation matrix
			 * @param[in] size_x : \f$n_x\f$, dimension of vector $\f$x\f$.
			 * - 0 OK
			 * - 1 problem
			 * @brief
			 * This function changes parameters.
			 */
			int set_params(const gsl_matrix * m,
						   unsigned int size_x);
		
		
			/**@fn void tkalman_equivalents :: free();
			 * @brief
			 * This function frees memory.
			 */
			void free();
		
			/**@fn int tkalman_equivalents :: alloc();
			 * @return
			 * - 0 ok
			 * @brief
			 * This function allocates object memory.
			 */
			int alloc();
		
			/**@fn void tkalman_equivalents :: create_views();
			 * @brief
			 * This function creates matrix views.
			 */
			void create_views();
		
		
			/**@fn void tkalman_equivalents :: initialize();
			 * @brief
			 * This function sets object attributes to 0.
			 */
			void initialize();
		//Données propres
			unsigned int _size_x;
			unsigned int _size_y;
			unsigned int _size_t;
			gsl_matrix * mat_tt;
			gsl_matrix * mat_tt_bis;
			gsl_matrix mat_xx;
			gsl_vector * vect_t;
			gsl_vector vect_x;
			gsl_permutation * perm_x;
			
		//Paramètres
			const gsl_matrix * m;
			gsl_matrix m_xx;
			gsl_matrix m_xy;
			
			//Paramètres
			gsl_matrix * m_inv;
			gsl_matrix m_inv_xx;
			gsl_matrix m_inv_xy;
	};



#endif
