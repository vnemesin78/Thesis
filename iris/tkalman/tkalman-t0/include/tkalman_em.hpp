/**@file tkalman_em.hpp
 * @author Valérian Némesin
 * @brief
 * Pairwise Kalman EM
 */

#ifndef _TKALMAN_EM42_HPP_
	#define _TKALMAN_EM42_HPP_
	
	#include <cmath>
	#include <gsl/gsl_matrix.h>
	#include <gsl/gsl_linalg.h>
	#include <gsl/gsl_blas.h>
	#include "gsl_triangle_matrix.hpp"
	#include "tkalman_sums.hpp"
	#include "tkalman_smoothing.hpp"
	#include "tkalman_prediction.hpp"
	#include "tkalman_filtering.hpp"
	#include "tkalman_equivalent.hpp"
	#include "tkalman_constants.hpp"
	#include "tkalman_argmax.hpp"
	

	/**@def tkalman_esperance_unref(moment, nb)
	 * @param moment : moment
	 * @param nb : nombre
	 * @brief
	 Cette macro désalloue un moment statistique vectoriel.
	 *
	**/
	#ifndef tkalman_expectation_unref

		#define tkalman_expectation_unref(moment, nb) if (moment) \
													{\
														for( unsigned int i = 0; i < (nb) ; ++ i) \
														{\
															if ( (moment)[i] )\
															{\
																gsl_vector_free((moment)[i]);\
															}\
														}\
														delete [] (moment);\
													}
	#endif

	#ifndef tkalman_covariance_unref

		/**@def tkalman_covariance_unref(moment, nb)
		 * @param moment : moment
		 * @param nb : nombre
		 * @brief
		 Cette macro désalloue un moment statistique matriciel.
		 *
		 */
		#define tkalman_covariance_unref(moment, nb) if (moment) \
													 {\
														for( unsigned int i = 0; i < (nb) ; ++ i) \
														{\
															if ( (moment)[i] )\
															{\
																gsl_matrix_free((moment)[i]);\
															}\
														}\
														delete [] (moment);\
													 }
	#endif
		
	#ifndef tkalman_expectation_ref_v2

		#define tkalman_expectation_ref_v2(moment, nb, size_x, size_t) if (! (moment) ) \
		{\
			(moment) = (gsl_vector **) malloc( sizeof(gsl_vector * ) * (nb) );\
			if (moment)\
			{\
				(moment)[0] = gsl_vector_calloc(size_t);\
				for (unsigned j = 1 ; j < (nb); ++ j)\
				{\
					(moment)[j] = gsl_vector_calloc(size_x);\
				}\
			}\
		}
	#endif

	#ifndef tkalman_covariance_ref_v2

		#define tkalman_covariance_ref_v2(moment, nb, size_1, size_2, size_1_0, size_2_0) if (! (moment) ) \
		{\
			(moment) = (gsl_matrix **) malloc( sizeof(gsl_matrix * ) * (nb) );\
			if (moment)\
			{\
				(moment)[0] = gsl_matrix_calloc(size_1_0, size_2_0);\
				for (unsigned j = 1 ; j < (nb) ; ++ j)\
				{\
					(moment)[j] = gsl_matrix_calloc(size_1, size_2);\
				}\
			}\
		}
	#endif

	/**@class tkalman_em
	 * @author 
     * Valérian Némesin
	 * @brief
	 * Pairwise Kalman EM algo.
	 * 
	 */
	class tkalman_em
	{
		public:
			/**@fn tkalman_em :: tkalman_em(	const gsl_vector * t_0,
											    const gsl_matrix * sqrt_q_0,
											    const gsl_matrix * f,
											    const gsl_matrix * sqrt_q,
											    unsigned int size_x,
											    unsigned int n,
											    unsigned int p);
			 * @param[in] t_0 : \f$ \hat{t}_0 \f$, initial state expectation
			 * @param[in] sqrt_q_0 : \f$ [Q_0]^{\frac{1}{2}} \f$, square root of initial state covariance matrix
			 * @param[in] f : \f$ F \f$, transition matrix
			 * @param[in] sqrt_q : \f$ [Q]^{\frac{1}{2}} \f$, square root of covariance matrix
			 * @param[in] n : number of samples
			 * @param[in] p : number of EM iterations
			 */
			tkalman_em(const gsl_vector * t_0,
					   const gsl_matrix * sqrt_q_0,
					   const gsl_matrix * f,
					   const gsl_matrix * sqrt_q,
					   unsigned int size_x,
					   unsigned int n,
					   unsigned int p);
			/**@fn int tkalman_em :: setup(	const gsl_vector * t_0,
											const gsl_matrix * sqrt_q_0,
											const gsl_matrix * f,
											const gsl_matrix * sqrt_q,
											unsigned int size_x,
											unsigned int n,
											unsigned int p);
			 * @param[in] t_0 : \f$ \hat{t}_0 \f$, initial state expectation
			 * @param[in] sqrt_q_0 : \f$ [Q_0]^{\frac{1}{2}} \f$, square root of initial state covariance matrix
			 * @param[in] f : \f$ F \f$, transition matrix
			 * @param[in] sqrt_q : \f$ [Q]^{\frac{1}{2}} \f$, square root of covariance matrix
			 * @param[in] n : number of samples
			 * @param[in] p : number of EM iterations
			 * @return
			 * - 0 OK
			 * - 1 else
			 * @brief
			 * Setup
			 */
			int setup(const gsl_vector * t_0,
					  const gsl_matrix * sqrt_q_0,
					  const gsl_matrix * f,
					  const gsl_matrix * sqrt_q,
					  unsigned int size_x,
					  unsigned int n,
					  unsigned int p);

			/**@fn bool tkalman_em :: operator ! () const;
			 * @return
			 * - 0 OK
			 * - 1 else
			 */
			virtual bool operator ! () const;
		
		
		
		
			/**@fn void tkalman_em :: filter( const gsl_vector * const * observations,
											 unsigned int n);
			 * @param[in] observations : observations
			 * @param[in] n : number of samples
			 * @brief
			 * filtered data.
			 * @warning
			 * Don't try \f$ n > n_0 \f$. 
			 */
			void filter(const gsl_vector * const * observations,
						unsigned int n);
		
		
		
			/**@fn void tkalman_em :: smooth( const gsl_vector * const * observations,
											  unsigned int n);
			 * @param[in] observations : observations
			 * @param[in] n : number of samples
			 * @brief
			 * smooth data.
			 * @warning
			 * Don't try \f$ n > n_0 \f$. 
			 */
			void smooth(const gsl_vector * const * observations,
						unsigned int n);
		
			//Accesseurs
			/**@fn inline const gsl_vector * tkalman_em :: t_0() const
			 * @return
			 * \f$ \hat{t}_0 \f$, initial state expectation
			 */
			inline const gsl_vector * t_0() const
			{
				return _t_0;
			}
			
			/**@fn inline const gsl_matrix * tkalman_em :: sqrt_q_0() const
			 * @return
			 * \f$ [Q_0]^{\frac{1}{2}} \f$, square root of initial state covariance matrix
			 */
			inline const gsl_matrix * sqrt_q_0() const
			{
				return _sqrt_q_0;
			}
			
			/**@fn void tkalman_em :: get_q_0(gsl_matrix * q_0) const
			 * @param[out] q_0: \f$Q_0\f$, initial state covariance matrix
			 */
			void get_q_0(gsl_matrix * q_0) const;
			
			/**@fn inline const gsl_matrix * tkalman_em :: f() const
			 * @return 
			 \f$ F \f$, transition matrix
			 */
			inline const gsl_matrix * f() const
			{
					return _f;
			}
			
			/**@fn inline const gsl_matrix * tkalman_em :: sqrt_q() const
			 * @return
			 * \f$ [Q]^{\frac{1}{2}} \f$, square root of covariance matrix
			 */
			inline const gsl_matrix * sqrt_q() const
			{
				return _sqrt_q;
			}
			 
			/**@fn void tkalman_em :: get_q(gsl_matrix * q) const
			 * @param[out] q: \f$ Q \f$, covariance matrix
			 */
			void get_q(gsl_matrix * q) const;
			
			
			/**@fn inline unsigned int tkalman_em :: size_x() const
			 * @return 
			 * \f$n_x\f$, dimension of x
			 */
			inline unsigned int size_x() const
			{
				return _size_x;
			}
			
			/**@fn inline unsigned int tkalman_em :: size_y() const
			 * @return 
			 * \f$n_y\f$, dimension of y
			 */
			inline unsigned int size_y() const
			{
				return _size_y;
			}
			
			/**@fn inline unsigned int tkalman_em :: size_t() const
			 * @return 
			 * \f$n_t\f$, dimension of t
			 */
			inline unsigned int size_t() const
			{
				return _size_t;
			}
			
			/**@fn const gsl_vector * const * tkalman_em :: x_p() const
			 * @return \f$ \hat{x}_{n|n-1}\f$, predicting state expectations
			 */
			inline const gsl_vector * const * x_p() const
			{
				return _x_p;
			}
			
			/**@fn const gsl_vector * const * tkalman_em :: x_f() const
			 * @return \f$ \hat{x}_{n|n}\f$, filtering state expectations
			 */
			inline const gsl_vector * const * x_f() const
			{
				return _x_f;
			}
			
			/**@fn const gsl_vector * const * tkalman_em :: x_s() const
			 * @return 
			 * \f$ \hat{x}_{n|N}\f$, smoothing state expectations
			 */
			inline const gsl_vector * const * x_s() const
			{
				return _x_s;
			}
			
			/**@fn inline const gsl_matrix * const * tkalman_em :: sqrt_p_p() const
			 * @return 
			 * \f$[P_{n|n-1}]^{\frac{1}{2}}\f$, square roots of predicting state covariance
			 */
			inline const gsl_matrix * const * sqrt_p_p() const
			{
				return _sqrt_p_p;
			}
			
			/**@fn inline const gsl_matrix * const * tkalman_em :: sqrt_p_s() const
			 * @return 
			 * \f$[P_{n|N}]^{\frac{1}{2}}\f$, square roots of smoothing state covariance
			 */
			inline const gsl_matrix * const * sqrt_p_s() const
			{
				return _sqrt_p_s;
			}
			
			/**@fn inline const gsl_matrix * const * tkalman_em :: sqrt_p_f() const
			 * @return 
			 * \f$[P_{n|n}]^{\frac{1}{2}}\f$, square roots of filtering state covariance
			 */
			inline const gsl_matrix * const * sqrt_p_f() const
			{
				return _sqrt_p_f;
			}
			//Nombre d'obs
			
			/**@fn inline unsigned int tkalman_em :: n() const
			 * @return 
			 * \f$N\f$, number of observations
			 */
			inline unsigned int n() const
			{
				return _n;
			}
			
			/**@fn inline unsigned int tkalman_em :: p() const
			 * @return 
			 * \f$p\f$, number of EM iterations
			 */
			inline unsigned int p() const
			{
				return _p;
			}
		
			/**@fn tkalman_em :: ~tkalman_em();
			 */
			~tkalman_em();
			
			
			/**@fn int tkalman_em :: set_params(const gsl_vector * t_0,
										   const gsl_matrix * sqrt_q_0,
										   const gsl_matrix * f,
										   const gsl_matrix * sqrt_q );
			 * @param[in] t_0 : \f$ \hat{t}_0 \f$, initial state expectation
			 * @param[in] sqrt_q_0 : \f$ [Q_0]^{\frac{1}{2}} \f$, square root of initial state covariance matrix
			 * @param[in] f : \f$ F \f$, transition matrix
			 * @param[in] sqrt_q : \f$ [Q]^{\frac{1}{2}} \f$, square root of covariance matrix
			 * @brief
			 * This function changes filter parameters.
			 */
			void set_params(const gsl_vector * t_0,
							const gsl_matrix * sqrt_q_0,
							const gsl_matrix * f,
							const gsl_matrix * sqrt_q);
							
			
		protected:
			/**@fn void tkalman_em :: initialize();
			 * @fn void tkalman_em :: initialize_params();
			 * @fn void tkalman_em :: initialize_objects();
			 * @fn void tkalman_em :: initialize_moments();
			 * @brief
			 * These methods set object attributes to 0.
			 */
			void initialize();
			void initialize_params();
			void initialize_objects();
			void initialize_moments();
		
			/**@fn int tkalman_em :: alloc();
			 * @fn int tkalman_em :: alloc_params();
			 * @fn int tkalman_em :: alloc_objects();
			 * @fn int tkalman_em :: alloc_moments();
			 * @return
			 * - 0 OK
			 * - 1 memory error
			 * @brief
			 * These methods allocate object memory.
			 */
			int alloc();
			int alloc_params();
			int alloc_objects();
			int alloc_moments();
			
			/**@fn void tkalman_em :: free();
			 * @fn void tkalman_em :: free_params();
			 * @fn void tkalman_em :: free_objects();
			 * @fn void tkalman_em :: free_moments();
			 * @brief
			 * These methods free object memory.
			 **/
			void free();
			void free_params();
			void free_objects();
			void free_moments();
			
			/**@fn bool tkalman_em :: check_params() const;
			 * @fn bool tkalman_em :: check_objects() const;
			 * @fn bool tkalman_em :: check_moments() const;
			 * @brief
			 * These methods check object attributes.
			 */
			bool check_params() const;
			bool check_objects() const;
			bool check_moments() const;
		
			/**@fn void tkalman_em :: create_views();
			 * @fn void tkalman_em :: create_param_views();
			 * @fn void tkalman_em :: create_moment_views();
			 * @brief
			 * These methods create different matrix views.
			 **/
			void create_views();
			void create_param_views();
			void create_moment_views();
			
			/**@fn void tkalman_em :: constraint()
			 * @brief
			 * This method constraints filter parameters. Likelihood is saved but an equivalent system where 
			 * \f$  =
			 * \begin{pmatrix}
			 * F^{x,x}	&	F^{x,y}	\\
			 * I		& 0
			 * \end{pmatrix}
			 * \f$ is found.
			 **/
			void constraint();
			

		
		
			//Objets
			tkalman_constants * constants;
			tkalman_prediction * prediction;
			tkalman_filtering * filtering;
			tkalman_smoothing * smoothing;
			tkalman_equivalents * equivalents;
			tkalman_sums * sums;
			tkalman_argmax * arg_max;
			
			//Paramètres
			gsl_vector * _t_0;
				gsl_vector x_0;
				gsl_vector y_m1;
			gsl_matrix * _sqrt_q_0;
			gsl_matrix * _f;
				gsl_matrix f_x_; // Vue sur la rangée x
				gsl_matrix f_y_; // Vue sur la rangée y
			gsl_matrix * _sqrt_q;
			unsigned int _n;
			unsigned int _p;
			unsigned int _size_x;
			unsigned int _size_y;
			unsigned int _size_t;
			unsigned int _i;
			//Constantes
			gsl_matrix * f2_x_,
					   * sqrt_q2_xx,
					   * q2_xy,
					   * sqrt_q_yy;
			//Equivalents
			gsl_matrix * _m;
				gsl_matrix m_x_; // Vue sur M
			
			//Moments
				gsl_vector ** _innovation;
				gsl_matrix ** _sqrt_s;
			
				gsl_vector ** _x_p;
				gsl_matrix ** _sqrt_p_p;
				gsl_vector ** _x_f;
		
				gsl_matrix ** _sqrt_p_f;
				gsl_vector ** _x_s;
				gsl_matrix ** _sqrt_p_s;
				gsl_matrix ** _c_s;
			//
				gsl_vector x_p_0;
				gsl_vector y_p_m1;
				gsl_vector x_f_0;
				gsl_vector y_f_m1;
				gsl_vector x_s_0;
				gsl_vector y_s_m1;
				
	};
	
	
	
	
	
	
	
#endif
