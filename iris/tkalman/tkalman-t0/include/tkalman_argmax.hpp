/**@file tkalman_argmax.hpp
 * @author Valérian Némesin
 * @brief
Pairwise Kalman maximisation
 *
 */
#ifndef _TKALMAN_ARGMAX2_HPP_
	#define _TKALMAN_ARGMAX2_HPP_
	#include <gsl/gsl_matrix.h>
	#include <gsl/gsl_linalg.h>
	#include <gsl/gsl_blas.h>
	#include <cmath>
	#include "gsl_triangle_matrix.hpp"

	/**@fn void tkalman_compute_argmax(gsl_vector * t_0,
									   gsl_matrix * sqrt_q_0,
									   gsl_matrix * f,
									   gsl_matrix * sqrt_q,
									   const gsl_matrix * sqrt_c_00,
									   const gsl_matrix * sqrt_c_01,
									   const gsl_matrix * sqrt_c_11,
									   unsigned int n,
									   const gsl_vector * t_s_0,
									   const gsl_matrix * sqrt_q_s_0,
									   gsl_matrix * mat_tt,
									   gsl_permutation * perm_t)
	 * @param[out] t_0 : \f$ \hat{t}_0 \f$, initial state expectation (NULL = no estimation)
	 * @param[out] sqrt_q_0 : \f$ [Q_0]^{\frac{1}{2}} \f$, square root of initial state covariance matrix (NULL = no estimation)
	 * @param[out] f : \f$ F \f$, transition matrix
	 * @param[out] sqrt_q : \f$ [Q]^{\frac{1}{2}} \f$, square root of covariance matrix
	 * @param[in] sqrt_c_00: \f$ [M_{sums}]^{\frac{1}{2}}_{0,0} \f$
	 * @param[in] sqrt_c_10 : \f$ [M_{sums}]^{\frac{1}{2}}_{1,0} \f$
	 * @param[in] sqrt_c_11 : \f$ [M_{sums}]^{\frac{1}{2}}_{1,1} \f$
	 * @param[in] n : number of samples
	 * @param[in] t_s_0 : \f$ \hat{t}_{0|N} \f$, initial smoothing state expectation (NULL = no estimation)
	 * @param[in] sqrt_q_s_0 : \f$ [Q_{0|N}]^{\frac{1}{2}} \f$, square root of initial smoothing state covariance matrix (NULL = no estimation)
	 * @param mat_tt : allocated \f$ ( n_t, n_t) \f$-matrix
	 * @param perm_t : allocated \f$ n_t-\f$permutation.
	 * @brief
	 This function computes the arguments of the maximum of auxiliary likelihood function:
	 - \f$ \hat{t}_0 = \hat{t}_{0|N} \f$
	 - \f$ Q_0 = Q_{0|N} \f$
	 - \f$ F = \tilde{C}_{1,0} [\tilde{C}_{0,0}]^{-1}\f$
	 - \f$ Q = \frac{1}{N+1} [\tilde{C}_{1,1} - \tilde{C}_{1,0} [\tilde{C}_{0,0}]^{-1} \tilde{C}_{0,1}]\f$
	 - \f$ M_{sums} = 
		\begin{pmatrix} 
			\tilde{C}_{0,0} 	& 	\tilde{C}_{0,1} \\ 
			\tilde{C}_{1,0}		& 	\tilde{C}_{1,1}
		\end{pmatrix}
		\f$
	 */
	void tkalman_compute_argmax(gsl_vector * t_0,
								gsl_matrix * sqrt_q_0,
								gsl_matrix * f,
								gsl_matrix * sqrt_q,
								const gsl_matrix * sqrt_c_00,
								const gsl_matrix * sqrt_c_01,
								const gsl_matrix * sqrt_c_11,
								unsigned int n,
								const gsl_vector * t_s_0,
								const gsl_matrix * sqrt_q_s_0,
								gsl_matrix * mat_tt,
								gsl_permutation * perm_t);


	/**@class
 * @author 
 * Valérian Némesin
	 * @brief
	 * Maximisation of auxiliary likelihood function.
	 **/
	class tkalman_argmax
	{
		public :
			/**@fn tkalman_argmax :: tkalman_argmax (const gsl_matrix * sqrt_c_00,
												     const gsl_matrix * sqrt_c_01,
												     const gsl_matrix * sqrt_c_11);
			 * @param[in] sqrt_c_00: \f$ [M_{sums}]^{\frac{1}{2}}_{0,0} \f$
			 * @param[in] sqrt_c_10 : \f$ [M_{sums}]^{\frac{1}{2}}_{1,0} \f$
			 * @param[in] sqrt_c_11 : \f$ [M_{sums}]^{\frac{1}{2}}_{1,1} \f$
			 * @brief
			 * Constructor
			 * - \f$ M_{sums} = 
				\begin{pmatrix} 
					\tilde{C}_{0,0} 	& 	\tilde{C}_{0,1} \\ 
					\tilde{C}_{1,0}		& 	\tilde{C}_{1,1}
				\end{pmatrix}
				\f$
			 * 
			 */
			tkalman_argmax(const gsl_matrix * sqrt_c_00,
						   const gsl_matrix * sqrt_c_01,
						   const gsl_matrix * sqrt_c_11);

			/**@fn int tkalman_argmax :: setup( const gsl_matrix * sqrt_c_00,
											    const gsl_matrix * sqrt_c_01,
											    const gsl_matrix * sqrt_c_11);
			 * @param[in] sqrt_c_00: \f$ [M_{sums}]^{\frac{1}{2}}_{0,0} \f$
			 * @param[in] sqrt_c_10 : \f$ [M_{sums}]^{\frac{1}{2}}_{1,0} \f$
			 * @param[in] sqrt_c_11 : \f$ [M_{sums}]^{\frac{1}{2}}_{1,1} \f$
			 * @return
			 * - 0 if success
			 * @brief
			 * Setup
			 * - \f$ M_{sums} = 
				\begin{pmatrix} 
					\tilde{C}_{0,0} 	& 	\tilde{C}_{0,1} \\ 
					\tilde{C}_{1,0}		& 	\tilde{C}_{1,1}
				\end{pmatrix}
				\f$
			**/
			virtual int setup(const gsl_matrix * sqrt_c_00,
							  const gsl_matrix * sqrt_c_01,
						      const gsl_matrix * sqrt_c_11);

			/**@fn inline void tkalman_argmax :: compute_argmax(	gsl_vector * t_0,
																	gsl_matrix * sqrt_q_0,
																	gsl_matrix * f,
																	gsl_matrix * sqrt_q,
																	unsigned int n,
																	const gsl_vector * t_s_0,
																	const gsl_matrix * sqrt_q_s_0);
			 * @param[out] t_0 : \f$ \hat{t}_0 \f$, initial state expectation (NULL = no estimation)
			 * @param[out] sqrt_q_0 : \f$ [Q_0]^{\frac{1}{2}} \f$, square root of initial state covariance matrix (NULL = no estimation)
			 * @param[out] f : \f$ F \f$, transition matrix
			 * @param[out] sqrt_q : \f$ [Q]^{\frac{1}{2}} \f$, square root of covariance matrix
			 * @param[in] n : number of samples
			 * @param[in] t_s_0 : \f$ \hat{t}_{0|N} \f$, initial smoothing state expectation (NULL = no estimation)
			 * @param[in] sqrt_q_s_0 : \f$ [Q_{0|N}]^{\frac{1}{2}} \f$, square root of initial smoothing state covariance matrix (NULL = no estimation)
			 * @brief
			 * This function computes the arguments of the maximum of auxiliary likelihood function:
			 * - \f$ \hat{t}_0 = \hat{t}_{0|N} \f$
			 * - \f$ Q_0 = Q_{0|N} \f$
			 * - \f$ F = \tilde{C}_{1,0} [\tilde{C}_{0,0}]^{-1}\f$
			 * - \f$ Q = \frac{1}{N+1} [\tilde{C}_{1,1} - \tilde{C}_{1,0} [\tilde{C}_{0,0}]^{-1} \tilde{C}_{0,1}]\f$
			 * - \f$ M_{sums} = 
				\begin{pmatrix} 
					\tilde{C}_{0,0} 	& 	\tilde{C}_{0,1} \\ 
					\tilde{C}_{1,0}		& 	\tilde{C}_{1,1}
				\end{pmatrix}
				\f$
			**/
			inline void compute_argmax(gsl_vector * t_0,
									   gsl_matrix * sqrt_q_0,
									   gsl_matrix * f,
									   gsl_matrix * sqrt_q,
									   unsigned int n,
									   const gsl_vector * t_s_0,
									   const gsl_matrix * sqrt_q_s_0)
			{
				tkalman_compute_argmax(t_0,
									   sqrt_q_0,
									   f,
									   sqrt_q,
									   _sqrt_c_00,
									   _sqrt_c_01,
									   _sqrt_c_11,
									   n,
									   t_s_0,
									   sqrt_q_s_0,
									   mat_tt,
									   perm_t);
			}

			/**@fn tkalman_argmax :: ~tkalman_argmax()
			 * @brief
			 * Destructor
			 */
			virtual ~tkalman_argmax();

			/**@fn bool tkalman_argmax :: operator!() const;
			 * @return
			 * - 0 no problem
			 * - 1 else
			 * @brief
			 * check object memory.
			 */
			virtual bool operator!() const;


		protected:
			/**@fn int tkalman_argmax :: set_params( 	const gsl_matrix * sqrt_c_00,
														const gsl_matrix * sqrt_c_01,
														const gsl_matrix * sqrt_c_11);
			 * @param[in] sqrt_c_00: \f$ [M_{sums}]^{\frac{1}{2}}_{0,0} \f$
			 * @param[in] sqrt_c_10 : \f$ [M_{sums}]^{\frac{1}{2}}_{1,0} \f$
			 * @param[in] sqrt_c_11 : \f$ [M_{sums}]^{\frac{1}{2}}_{1,1} \f$
			 * @brief
			 * This function changes the sums.
			 */
			int set_params(const gsl_matrix * sqrt_c_00,
						   const gsl_matrix * sqrt_c_01,
						   const gsl_matrix * sqrt_c_11);

			/**@fn void tkalman_argmax :: free();
			 * @brief
			 * This function frees memory.
			 */
			void free();

			/**@fn int tkalman_argmax :: alloc();
			 * @return
			 * - 0 no problem
			 * @brief
			 * Memory allocation.
			 */
			int alloc();

			/**@fn void tkalman_argmax :: create_views();
			 * @brief
			 * This function creates matrix views.
			 */
			inline void create_views()
			{}

			/**@fn void tkalman_argmax :: initialize();
			 * @brief
			 * This function sets the object attributes to 0.
			 */
			void initialize();


		//Données propres
			unsigned int _size_t;
			gsl_matrix * mat_tt;
			gsl_permutation * perm_t;
		//Paramètres
			const gsl_matrix * _sqrt_c_00;
			const gsl_matrix * _sqrt_c_01;
			const gsl_matrix * _sqrt_c_11;
	};

#endif
