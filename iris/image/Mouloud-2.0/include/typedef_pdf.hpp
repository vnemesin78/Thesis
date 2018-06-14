/**@file typedef_pdf.hpp
 * 
 */
#ifndef _TYPEDEF_PDF_HPP_
#define _TYPEDEF_PDF_HPP_
	
	/**@typedef 
	 * @brief
	 * Typedef des densités de proba
	 */
	typedef double (*pdf_type) ( const double & x,
									const void * params );
									
	typedef double (*pdf_field_type ) (	const double & _r,
											const double & r,
											const double & r_,
											const void * params );
									
									
									
	typedef int (*estimation_type)  (  void * params,
										  const double * data,
										  const unsigned char * mask,
										  unsigned int width,
										  unsigned int height,
										  unsigned int width_step,
										  unsigned int mask_width_step,
										  unsigned char mask_value,
										  const void * other_params );
#endif
