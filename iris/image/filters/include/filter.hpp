#ifndef _FILTER_HPP_
	#define _FILTER_HPP_
	#include <fftw3.h>
	
	/**@fn void convolve (	fftw_complex * signal_out,
							const fftw_complex * signal_in,
							const fftw_complex * filter_fft,
							unsigned int width,
							fftw_complex * tmp )
	 * @param signal_out : signal filtré
	 * @param signal_in : signal original
	 * @param filter_fft : réponse impulsonielle du filtre
	 * @param width : taille des signaux
	 * @brief
	 * Convolution 
	 **/
	void convolve (	fftw_complex * signal_out,
					fftw_complex * signal_in,
					const fftw_complex * filter_fft,
					unsigned int width,
					fftw_complex * tmp );
	
	/**@fn void filter_1d( 	double * signal_out,
							const double * signal_in,
							unsigned int size,
							const double * kernel,
							unsigned int kernel_size )
	 * @param signal_out : signal filtré
	 * @param signal_in : signal d'origine
	 * @param size : taille du signl
	 * @param kernel : noyau
	 * @param kernel_size : taille du noyau
	 * 
	 */
	void filter_1d( double * signal_out,
					const double * signal_in,
					unsigned int size,
					const double * kernel,
					unsigned int kernel_size );
					
	void convolve (	const fftw_complex * filter_fft,
					unsigned int width,
					fftw_complex * tmp,
					fftw_plan & plan_in,
					fftw_plan & plan_out );
					
	void convolve_build_fftw_plans (	fftw_plan & plan_in,
										fftw_plan & plan_out,
										fftw_complex * signal_out,
										fftw_complex * signal_in,
										unsigned int width,
										fftw_complex * tmp );
#endif
