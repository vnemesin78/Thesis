/**@file gibbs_pdf.hpp
 * 
 * 
 */
#ifndef _GIBBS_PDF_HPP_
#define _GIBBS_PDF_HPP_
	
	/**@struct
	 * @brief
	 * Paramètres d'un champs de gibbs
	 * 
	 */
	struct gibbs_params
	{
		double a, delta;
	};


	
	
	/**@fn
	 * @param x : r_i
	 * @param params:
	 * -arg[0] : \f$ \delta \f$
	 * -arg[1] : \f$ r_{i-1} \f$
	 * -arg[2] : \f$ r_{i+1} \f$
	 * @brief
	 * Computes:
	 * \f$ A \times \exp( \delta ( (r_i - r_{i-1} )^2 + (r_{i+1} - r_{i} )^2 ) )\f$
	 */
	double gibbs_pdf( 	const double & _r,
						const double & r,
						const double & r_,
						const void * params );
	
	
	
	/**@fn
	 * @param x : variable
	 * @param params:
	 * -arg[0] : \f$ \delta \f$
	 * -arg[1] : \f$ r_{i-1} \f$
	 * -arg[2] : \f$ r_{i+1} \f$
	 * @brief
	 * Computes:
	 * \f$ \delta ( (r_i - r_{i-1} )^2 + (r_{i+1} - r_{i} )^2 )\f$
	 */
	double log_gibbs_pdf( 	const double & _r,
							const double & r,
							const double & r_,
							const void * params );
	
	
	
	
	
	
	
#endif
