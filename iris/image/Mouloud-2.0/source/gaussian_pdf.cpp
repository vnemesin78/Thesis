#include "gaussian_pdf.hpp"
#include <cmath>
#include <iostream>
using namespace std;
double gaussian_pdf ( 	const double & x,
						const void * params )
{
	gaussian_params * data = (gaussian_params*) params;
	double dx = (x - data->mean ) / data->sigma;
	return ( exp ( - dx * dx ) / ( sqrt( 2 * M_PI) * data->sigma ) );
}


double log_gaussian_pdf ( 	const double & x,
							const void * params )
{
	gaussian_params * data = (gaussian_params*) params;
	double dx = (x - data->mean ) / data->sigma;
	return ( -( log ( sqrt( 2 * M_PI) * data->sigma  ) + dx * dx ) );
}

int gaussian_estimation ( void * params,
						  const double * image,
						  const unsigned char * mask,
						  unsigned int width,
						  unsigned int height,
						  unsigned int width_step,
						  unsigned int mask_width_step,
						  unsigned char mask_value,
						  const void * other_params )
{
	gaussian_params * data = (gaussian_params*) params;
	//Calcul des moments
	double m1 = 0,
			m2 = 0;
	unsigned int nb_pixels = 0;
	
	for ( unsigned int i = 0; i < height; ++ i )
	{
		for ( unsigned int j = 0; j < width; ++ j )
		{
			if ( mask[ i * mask_width_step + j ] == mask_value )
			{
				double v = image[ i * width_step + j ];
				nb_pixels ++;
				m1 += v;
				m2 += v * v;
			}
		}
	}
	if ( nb_pixels == 0 )
		return 1;
	
	data->mean = m1 / nb_pixels;
	data->sigma = sqrt( m2 / nb_pixels - data->mean * data->mean );
	return 0;
}



