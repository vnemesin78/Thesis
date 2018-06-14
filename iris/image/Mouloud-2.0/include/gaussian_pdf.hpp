/**@file gaussian_pdf.hpp
 * 
 * 
 */
#ifndef _GAUSSIAN_PDF_HPP_
#define _GAUSSIAN_PDF_HPP_

	/**@struct
	 * @brief
	 * Paramètres d'une loi gaussienne
	 * 
	 */
	struct gaussian_params
	{
		double mean,
				sigma;
	};


	/**@fn
	 * @param x : variable
	 * @param params : law parameters
	 * - arg[0] -> mean
	 * - arg[1] -> variance
	 * @brief
	 * Gaussian pdf
	 */
	double gaussian_pdf ( const double & x,
						   const void * params );


	/**@fn
	 * @param x : variable
	 * @param params : law parameters
	 * - arg[0] -> mean
	 * - arg[1] -> variance
	 * @brief
	 * log-Gaussian pdf
	 */
	double log_gaussian_pdf ( const double & x,
							   const void * params );

	/**@fn
	 * @brief
	 * Updating of gaussian parameters.
	 * 
	 */
	int gaussian_estimation ( void * params,
							  const double * data,
							  const unsigned char * mask,
							  unsigned int width,
							  unsigned int height,
							  unsigned int width_step,
							  unsigned int mask_width_step,
							  unsigned char mask_value,
							  const void * other_params );

#endif
