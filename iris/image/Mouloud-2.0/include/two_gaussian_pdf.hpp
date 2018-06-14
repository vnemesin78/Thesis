/**@file two_gaussian_pdf.hpp
 * 
 * 
 */
#ifndef _TWO_GAUSSIAN_PDF_HPP_
#define _TWO_GAUSSIAN_PDF_HPP_
	#include "gaussian_pdf.hpp"

	/**@struct
	 * @brief
	 * Paramètres d'une loi gaussienne
	 * 
	 */
	struct two_gaussian_params
	{
		double weight_1,
				weight_2;
		gaussian_params law_1,
						  law_2;
	};

	struct two_gaussian_estimation_params
	{
		unsigned int nb_iterations_algo_em;
	};



	/**@fn
	 * @param x : variable
	 * @param params : law parameters
	 * - arg[0] -> weight 1
	 * - arg[1] -> mean 1
	 * - arg[2] -> variance 1
	 * - arg[3] -> weight 2
	 * - arg[4] -> mean 2
	 * - arg[5] -> variance 2
	 * @brief
	 * log-Gaussian pdf
	 */
	double two_gaussian_pdf(  const double & x,
							   const void * params );
	/**@fn
	 * @param x : variable
	 * @param params : law parameters
	 * - arg[0] -> weight 1
	 * - arg[1] -> mean 1
	 * - arg[2] -> variance 1
	 * - arg[3] -> weight 2
	 * - arg[4] -> mean 2
	 * - arg[5] -> variance 2
	 * @brief
	 * log-Gaussian pdf
	 */
	double log_two_gaussian_pdf(  	const double & x,
									const void * params );

	/**@fn
	 * @brief
	 * Estimation des paramètres de 2 lois gaussiennes par algo em.
	 */
	int two_gaussian_estimation ( void * params,
								  const double * data,
								  const unsigned char * mask,
								  unsigned int width,
								  unsigned int height,
								  unsigned int width_step,
								  unsigned int mask_width_step,
								  unsigned char mask_value,
								  const void * other_params );

#endif
