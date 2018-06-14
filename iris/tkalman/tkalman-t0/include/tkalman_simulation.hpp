/**@file tkalman_simulation.hpp
 * @author Valérian Némesin
 * @brief
 * Pairwise Kalman simulation
 */
#ifndef _TKALMAN_SIMULATION_2_HPP_
#define _TKALMAN_SIMULATION_2_HPP_
#include "tkalman_em.hpp"
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include "tkalman_experiment.hpp"
/**@class tkalman_simulation_2
 * @author 
 * Valérian Némesin
 * @brief
 * 
 * Pairwise Kalman simulation
 * 
 **/
class tkalman_simulation_2
{
	public:
		/**@fn tkalman_simulation_2 :: tkalman_simulation_2(	const gsl_vector * t0,
																const gsl_matrix * sqrt_q0,
																const gsl_matrix * f,
																const gsl_matrix * sqrt_q,
																gsl_rng * random_number_generator,
																unsigned int size_x,
																unsigned int nb_observations);
		 * @param[in] t0 : \f$ \hat{t}_0 \f$, initial state expectation
		 * @param[in] sqrt_q0 : \f$ [Q_0]^{\frac{1}{2}} \f$, square root of initial state covariance matrix
		 * @param[in] f : \f$ F \f$, transition matrix
		 * @param[in] sqrt_q : \f$ [Q]^{\frac{1}{2}} \f$, square root of covariance matrix
		 * @param random_number_generator : random number generator
		 * @param[in] size_x : dimension of vector \f$x\f$.
		 * @param[in] nb_observations : number of observations
		 */
		tkalman_simulation_2(const gsl_vector * t0,
							 const gsl_matrix * sqrt_q0,
							 const gsl_matrix * f,
							 const gsl_matrix * sqrt_q,
							 gsl_rng * random_number_generator,
							 unsigned int size_x,
							 unsigned int nb_observations);
							 
	
		/**@fn int tkalman_simulation_2 :: setup(	const gsl_vector * t0,
													const gsl_matrix * sqrt_q0,
													const gsl_matrix * f,
													const gsl_matrix * sqrt_q,
													gsl_rng * random_number_generator,
													unsigned int size_x,
													unsigned int nb_observations);
		 * @param[in] t0 : \f$ \hat{t}_0 \f$, initial state expectation
		 * @param[in] sqrt_q0 : \f$ [Q_0]^{\frac{1}{2}} \f$, square root of initial state covariance matrix
		 * @param[in] f : \f$ F \f$, transition matrix
		 * @param[in] sqrt_q : \f$ [Q]^{\frac{1}{2}} \f$, square root of covariance matrix
		 * @param random_number_generator : random number generator
		 * @param[in] size_x : dimension of vector \f$x\f$.
		 * @param[in] nb_observations : number of observations
		 * @return
		 * - 0 OK
		 * - 1 problem
		 * @brief
		 * Setup
		 */
		int setup(const gsl_vector * t0,
				  const gsl_matrix * sqrt_q0,
				  const gsl_matrix * f,
				  const gsl_matrix * sqrt_q,
				  gsl_rng * random_number_generator,
				  unsigned int size_x,
				  unsigned int nb_observations);
				
				
		/**@fn void tkalman_simulation_2 :: do_simulation();	
		 * @brief
		 * Do simulations.
		 **/  
		void do_simulation();		
		  
		/**@fn void tkalman_simulation_2 :: get_instant_mean_square_error(gsl_vector ** i_mse,
																		  const gsl_vector * x_est,
																		  unsigned int n,
																		  double scale);
		 * @param i_mse : instant mean square error or (i_mse = scale x i_mse(0) + instant mean square error)
		 * @param x_est : estimated states
		 * @param n : number of estimated states
		 * @param scale :
		 * @brief
		 * This function computes the instant mean square error of data.
		 */
		void get_instant_mean_square_error(gsl_vector ** i_mse,
										   const gsl_vector * const * x_est,
										   unsigned int n,
										   double scale = 0);
		
		/**@fn inline const gsl_vector * const * tkalman_simulation_2 :: x() const
		 * @return
		 * hidden states
		 **/
		inline const gsl_vector * const * x() const
		{
			return _x;
		}
		
		/**@fn inline const gsl_vector * const * tkalman_simulation_2 :: y() const
		 * @return
		 * observations
		 **/
		inline const gsl_vector * const * y() const
		{
			return _y;
		}
		/**@fn inline unsigned int tkalman_simulation_2 :: n() const
		 * @return
		 * number of observations
		 **/
		inline unsigned int n() const
		{
			return _n;
		}
		/**@fn inline unsigned int tkalman_simulation_2 :: size_x() const
		 * @return
		 * dimension of vector \f$x\f$.
		 **/
		inline unsigned int size_x() const
		{
			return _size_x;
		}
		/**@fn inline unsigned int tkalman_simulation_2 :: size_y() const
		 * @return
		 * dimension of vector \f$y\f$.
		 **/
		inline unsigned int size_y() const
		{
			return _size_y;
		}
		/**@fn inline unsigned int tkalman_simulation_2 :: size_t() const
		 * @return
		 * dimension of vector \f$t\f$.
		 **/
		inline unsigned int size_t() const
		{
			return _size_t;
		}

		
		/**@fn tkalman_simulation_2 :: ~tkalman_simulation_2();
		 * 
		 **/
		~tkalman_simulation_2();
	
	
		/**@fn bool tkalman_simulation_2 :: operator!() const;
		 * @return
		 * - 0 OK
		 * - 1 problem
		 * @brief
		 * This method checks object memory.
		 **/
		bool operator!() const;
	
	protected:
		
		void initialize();
		void initialize_moments();
		void initialize_tmp();
		void initialize_params();
		
		void free();
		void free_tmp();
		void free_moments();
		
		int alloc();
		int alloc_tmp();
		int alloc_moments();
		
		int check_params() const;
		int check_tmp() const;
		int check_moments() const;
		
		void create_tmp_view();
		
		
		//Moments
		gsl_vector ** _x;
		gsl_vector ** _y;
		//Tmp
		gsl_vector * vect_t_1;
				gsl_vector vect_t_1_view_x;
				gsl_vector vect_t_1_view_y;
		gsl_vector * vect_t_2;
		
		//Params
		unsigned int _size_x;
		unsigned int _size_y;
		unsigned int _size_t;
		
		const gsl_vector * _t0;
		const gsl_matrix * _sqrt_q0;
		const gsl_matrix * _f;
		const gsl_matrix * _sqrt_q;
		gsl_rng * rng;
		unsigned int _n;
};


/**@fn void do_tkalman_simulation_2(gsl_vector ** x,
									 gsl_vector ** y,
									 const gsl_vector * t0,
									 const gsl_matrix * sqrt_q0,
									 const gsl_matrix * f,
									 const gsl_matrix * sqrt_q,
									 const unsigned int nb_observations,
									 gsl_vector * vect_t_1,
									 gsl_vector * vect_t_1_view_x,
									 gsl_vector * vect_t_1_view_y,
									 gsl_vector * vect_t_2,
									 gsl_rng * r)
 * @param[out] x : hidden states ( n + 1 elements (0 to n) )
 * @param[out] y : observations ( n elements (0 to n - 1) )
 * @param[in] t0 : \f$ \hat{t}_0 \f$, initial state expectation
 * @param[in] sqrt_q0 : \f$ [Q_0]^{\frac{1}{2}} \f$, square root of initial state covariance matrix
 * @param[in] f : \f$ F \f$, transition matrix
 * @param[in] sqrt_q : \f$ [Q]^{\frac{1}{2}} \f$, square root of covariance matrix
 * @param[in] nb_observations : number of observations
 * @param vect_t_1 : \f$V\f$, allocated \f$n_t\f$-vector
 * @param vect_t_1_view_x : vector view on \f$V\f$ starting at 0, ending at n_x - 1.
 * @param vect_t_1_view_y : vector view on \f$V\f$ starting at n_x, ending at n_t
 * @param vect_t_2 : allocated \f$n_t\f$-vector
 * @param r :random number generator
 * @brief
 * This function simulates pairwise Kalman data.
**/
void do_tkalman_simulation_2(gsl_vector ** x,
						     gsl_vector ** y,
						     const gsl_vector * t0,
						     const gsl_matrix * sqrt_q0,
						     const gsl_matrix * f,
						     const gsl_matrix * sqrt_q,
						     const unsigned int nb_observations,
						     gsl_vector * vect_t_1,
						     gsl_vector * vect_t_1_view_x,
						     gsl_vector * vect_t_1_view_y,
						     gsl_vector * vect_t_2,
						     gsl_rng * r);
						     
/**@fn void do_no_tkalman_simulation_2(gsl_vector ** y,
									   const gsl_vector * const * x,
									   const gsl_matrix * sqrt_q_yy,
									   const unsigned int nb_observations,
									   gsl_vector * vect,
									   gsl_rng * r);

 * @param[out] y : observations ( n elements (0 to n - 1) )
 * @param[in] x : hidden states ( n + 1 elements (0 to n) )
 * @param[in]  sqrt_q_yy :
  \f$[Q^{y,y}]^{\frac{1}{2}}\f$, square root of the measurement covariance noise.
 * @param[in] nb_observations : number of observations
 * @param vect : \f$V\f$, allocated \f$n_x\f$-vector
 * @param r :random number generator
 * @brief
 * This function simulates data with a Gaussian noise.
 */
void do_no_tkalman_simulation_2(gsl_vector ** y,
								const gsl_vector * const * x,
								const gsl_matrix * sqrt_q_yy,
								const unsigned int nb_observations,
								gsl_vector * vect,
						        gsl_rng * r);

#endif
