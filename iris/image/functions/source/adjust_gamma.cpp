#include "adjust_gamma.hpp"
#include <cmath>
int adjust_gamma (	double * img_data,
					unsigned int width,
					unsigned int height,
					unsigned int width_step,
					double g )
{
	double 	min,
			max;
	double inv_g;
	//Gamma value must be > 0.
	if ( g <= 0 )
		return 1;
	
	inv_g = 1.0 / g;
	
	//Recherche du min
	min = img_data[0];
	for ( unsigned int i = 0; i < height; ++ i )
	{
		for ( unsigned int j = 0; j < width; ++ j )
		{
			unsigned int n_pixel = i * width_step + j;
			if ( img_data[n_pixel] < min )
				min = img_data[n_pixel];
		}
		
	}
	
	//Recherche du max
	max = img_data[0];
	for ( unsigned int i = 0; i < height; ++ i )
	{
		for ( unsigned int j = 0; j < width; ++ j )
		{
			unsigned int n_pixel = i * width_step + j;
			if ( img_data[n_pixel] > max )
				max = img_data[n_pixel];
		}
	}
	
	//Image uniforme
	if ( max == min )
		return 1;
	
	//Correction du gamma
	for ( unsigned int i = 0; i < height; ++ i )
	{
		for ( unsigned int j = 0; j < width; ++ j )
		{
			unsigned int n_pixel = i * width_step + j;
			img_data[n_pixel] = pow( ( img_data[n_pixel] - ( (double) min ) ) / ( max - ( (double) min ) ) , inv_g);
		}
	}
	return 0;
}

int adjust_gamma_2 ( 	double * img_data,
						unsigned int width,
						unsigned int height,
						unsigned int width_step,
						double g )
{
	double inv_g;
	//Gamma value must be > 0.
	if ( g <= 0 )
		return 1;
	
	inv_g = 1.0 / g;
	
	for ( unsigned int i = 0; i < height; ++ i )
	{
		for ( unsigned int j = 0; j < width; ++ j )
		{
			unsigned int n_pixel = i * width_step + j;
			img_data[n_pixel] =  pow( img_data[n_pixel] / 255.0 , inv_g);
		}
	}
	return 0;
	
}
