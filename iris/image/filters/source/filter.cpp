#include "filter.hpp"
#include <iostream>
using namespace std;





void filter_1d( double * signal_out,
				const double * signal_in,
				unsigned int size,
				const double * kernel,
				unsigned int kernel_size )
{
	unsigned int _h_k_size,
				 _r_k,
				 i,
				 j,
				 beg,
				 end;
	
	_h_k_size = (kernel_size - 1) / 2;
	_r_k = size - _h_k_size;
	//Filtrage à gauche
	for ( i = 0; i < size; ++ i )
	{
		signal_out[i] = 0;
		if ( i < _h_k_size )
			beg = _h_k_size - i;
		else
			beg = 0;
			
		if ( i >= _r_k )
		{
			end = size - i;
		}
		else
			end = kernel_size;
		for ( j = beg; j < end; ++ j)
		{
			signal_out[i] += signal_in[ j + i - _h_k_size ] * kernel[j];
		}
	}
}


void convolve (	fftw_complex * signal_out,
				fftw_complex * signal_in,
				const fftw_complex * filter_fft,
				unsigned int width,
				fftw_complex * tmp )
{
	fftw_plan p;
	unsigned int i;
	//FFTW de l'image
	p = fftw_plan_dft_1d(	width, 
							signal_in, 
							tmp, 
							FFTW_FORWARD, 
							FFTW_ESTIMATE);
	fftw_execute(p);
	fftw_destroy_plan(p);
	
	for ( i = 0; i < width; ++ i)
	{

		
		
		fftw_complex c;
		c[0] = tmp[i][0];
		c[1] = tmp[i][1];

		tmp[i][0] = c[0] * filter_fft[i][0] - c[1] * filter_fft[i][1];
		tmp[i][1] = c[1] * filter_fft[i][0] + c[0] * filter_fft[i][1];
	}
	
	//FFTW inverse de la réponse impu de l'image filtrée
	p = fftw_plan_dft_1d(	width, 
							tmp, 
							signal_out, 
							FFTW_BACKWARD, 
							FFTW_ESTIMATE);
	fftw_execute(p);
	fftw_destroy_plan(p);
}

void convolve (	const fftw_complex * filter_fft,
				unsigned int width,
				fftw_complex * tmp,
				fftw_plan & plan_in,
				fftw_plan & plan_out )
{
	unsigned int i;
	fftw_execute(plan_in);
	for ( i = 0; i < width; ++ i)
	{

		
		
		fftw_complex c;
		c[0] = tmp[i][0];
		c[1] = tmp[i][1];

		tmp[i][0] = c[0] * filter_fft[i][0] - c[1] * filter_fft[i][1];
		tmp[i][1] = c[1] * filter_fft[i][0] + c[0] * filter_fft[i][1];
	}
	
	//FFTW inverse de la réponse impu de l'image filtrée
	fftw_execute(plan_out);
	
}

void convolve_build_fftw_plans (	fftw_plan & plan_in,
									fftw_plan & plan_out,
									fftw_complex * signal_out,
									fftw_complex * signal_in,
									unsigned int width,
									fftw_complex * tmp )
{
	plan_in = fftw_plan_dft_1d(	width, 
								signal_in, 
								tmp, 
								FFTW_FORWARD, 
								FFTW_ESTIMATE);
	
	
	plan_out = fftw_plan_dft_1d(	width, 
									tmp, 
									signal_out, 
									FFTW_BACKWARD, 
									FFTW_ESTIMATE);
}
