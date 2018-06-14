#include "tkalman_experiment.hpp"

/**@fn void tkalman_get_instant_mean_square_error(gsl_vector ** i_mse,
												  const gsl_vector * const * x,
												  const gsl_vector * const * x_est,
												  const unsigned int n,
												  double scale,
												  gsl_vector * vect_x)
 * @param i_mse : instant mean square errors or (i_mse = scale x i_mse(0) + instant mean square error)
 * @param x : hidden states
 * @param x_est : estimated states
 * @param n : number of estimated states
 * @param scale : cf i_mse
 * @param vect_x : allowed vector of x size.
 * @brief
 * This function computes the instant mean square error of an sample of data.
 */
void tkalman_get_instant_mean_square_error(gsl_vector ** i_mse,
										   const gsl_vector * const * x,
										   const gsl_vector * const * x_est,
										   const unsigned int n,
										   double scale,
										   gsl_vector * vect_x)
{
	for (unsigned int i = 0 ; i < n ; i++)
	{
		gsl_vector_const_view view = gsl_vector_const_subvector(x[i], 0, vect_x->size);
		gsl_vector_const_view view2 = gsl_vector_const_subvector(x_est[i], 0, vect_x->size);
		gsl_vector_memcpy(vect_x, &(view2.vector));
		gsl_vector_sub(vect_x, &(view.vector));
		gsl_vector_mul(vect_x, vect_x);
		gsl_vector_scale(i_mse[i], scale);
		gsl_vector_add(i_mse[i], vect_x);
	}
}

