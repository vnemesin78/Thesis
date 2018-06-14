#include "gabor_filter.hpp"
#include <cmath>
#include <cstring>
void create_gabor_filter (	fftw_complex * gabor_filter,
							unsigned int width,
							double wavelenght,
							double sigma )
{
	unsigned int i,
				 width_s_2;
	double _2_log2_sigma_inv;
	double tmp;
	double radius;
	//Mise à zéro de la sortie
	memset (	gabor_filter,
				0,
				sizeof ( fftw_complex ) * width );
	
	
	//Calcul des constantes
	_2_log2_sigma_inv = log( sigma );
	_2_log2_sigma_inv = _2_log2_sigma_inv * _2_log2_sigma_inv * 2;
	_2_log2_sigma_inv = 1 / _2_log2_sigma_inv;
	
	width_s_2 = width / 2;
	
	
	//Calcul du filtre
	gabor_filter[0][0] = 0;  
	
	for ( i = 1; i < width_s_2 ; ++ i )
	{
		radius = ( 0.5 * i ) / width_s_2;
		tmp = log ( wavelenght * radius );
		tmp = tmp * tmp;
		gabor_filter[i][0] = exp( - tmp * _2_log2_sigma_inv );  
		
	}
}

