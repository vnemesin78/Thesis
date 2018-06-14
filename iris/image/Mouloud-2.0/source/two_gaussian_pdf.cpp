#include "two_gaussian_pdf.hpp"
#include "gaussian_pdf.hpp"
#include <cmath>
#include <cstring>
#include <cstdio>
#include <iostream>
using namespace std;

double two_gaussian_pdf(  	const double & x,
							const void * params )
{
	two_gaussian_params * data = (two_gaussian_params*) params;
	
	return ( data->weight_1 * gaussian_pdf( x, (void*) & ( data->law_1 ) )
		+  data->weight_2 * gaussian_pdf( x, (void*) & ( data->law_2 ) ) );
}

double log_two_gaussian_pdf(  	const double & x,
								const void * params )
{
	return log( two_gaussian_pdf( x, params ) );
}

int two_gaussian_estimation ( void * params,
							  const double * image,
							  const unsigned char * mask,
							  unsigned int width,
							  unsigned int height,
							  unsigned int width_step,
							  unsigned int mask_width_step,
							  unsigned char mask_value,
							  const void * other_params )
{
	two_gaussian_params * data = (two_gaussian_params*) params;
	unsigned int nb_iter = ( ( two_gaussian_estimation_params *) other_params )->nb_iterations_algo_em;
	
	unsigned char * mask_bis = new unsigned char[ width * height ];
	for ( unsigned int i = 0; i < nb_iter; ++ i )
	{
		memset( mask_bis,
				0,
				sizeof( char ) * width * height );
		double nb_pixels_1 = 0,
				nb_pixels_2 = 0;
		
		//Classification
		for ( unsigned int l = 0; l < height; ++ l )
		{
			for ( unsigned int j = 0; j < width; ++j )
			{
				if ( mask[ l * mask_width_step + j ] == mask_value )
				{
					double v = image[ l * width_step + j ];


					if ( data->weight_1 * gaussian_pdf( v, (void*) & ( data->law_1 ) ) >
						 data->weight_2 * gaussian_pdf( v, (void*) & ( data->law_2 ) ) )
					{
						mask_bis[l * width + j] = 1;
						nb_pixels_1 ++;
					}
					else
					{
						mask_bis[l * width + j] = 2;
						nb_pixels_2 ++;
					}
				}
				
				
			}
		}
		
		
		//Estimation des lois
		data->weight_1 = nb_pixels_1 / ( nb_pixels_1 + nb_pixels_2 );
		data->weight_2 = 1.0 - data->weight_1;
		
		int q = 0;
		if ( gaussian_estimation (  (void*) & ( data->law_1 ),
									image,
									mask_bis,
									width,
									height,
									width_step,
									width,
									1,
									NULL ) )
			q = 1;
		if ( gaussian_estimation (  (void*) & ( data->law_2 ),
									image,
									mask_bis,
									width,
									height,
									width_step,
									width,
									2,
									NULL ) )
			q ++;
		

		
		if ( q == 2 )
		{
			delete[] mask_bis;
			return 1;
		}
		else if ( q == 1 )
			break;
	}
	delete[] mask_bis;
	return 0;
}




