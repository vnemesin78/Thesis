/**@file tkalman_experiment.hpp
 * @brief
 * Statistical moments.
 * @author 
 * Valérian Némesin
 */
#ifndef _TKALMAN_EXPERIMENT_HPP_
#define _TKALMAN_EXPERIMENT_HPP_
#include "tkalman_em.hpp"
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
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
										   gsl_vector * vect_x);
	
#endif
