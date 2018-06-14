#include "chistogram.hpp"
#include <cstring>
#include <cmath>
#include <algorithm>
#include <iostream>

#include "filter.hpp" //A améliorer

using namespace std;
c_histogram :: c_histogram (	unsigned int nb_bins )
{
	initialize();
	setup ( nb_bins );
}

void c_histogram :: setup (	unsigned int nb_bins )
{
	if ( nb_bins == 0 )
		nb_bins = DEFAULT_NB_BINS;
	free();
	initialize();
	
	_nb_bins = nb_bins;
	
	_data = new double[ DEFAULT_NB_BINS ];
	_s_data = new double[ DEFAULT_NB_BINS ];
	_d_data = new double[ DEFAULT_NB_BINS];
	_d_d_data = new double[ DEFAULT_NB_BINS];
	
}

template <class type> void c_histogram :: compute_histogram (	const type * data,
																unsigned int width,
																unsigned int height,
																unsigned int width_step)
			
{
	unsigned int n_pixel;
	if ( width_step == 0)
		width_step = width;
		
	memset ( _data, 0, sizeof( double ) * _nb_bins ); // 0
	_nb_pixels = 0; // Nb pixels
	
	//Calcul des valeurs qui bornent l'histogramme
	for ( unsigned int i = 0; i < height; ++ i )
	{
		for ( unsigned int j = 0; j < width; ++ j )
		{
			n_pixel = i * width_step + j;
			if (  data[ n_pixel ] < _nb_bins )
			{

				_nb_pixels ++;
				_data[ (unsigned char) data[ n_pixel ] ]++;
			}
		}
	}
	

}

template <class type> void  c_histogram :: compute(	const type * data,
														const unsigned char * mask_data,
														unsigned int width,
														unsigned int height,
														unsigned int width_step,
														unsigned int mask_width_step,
														const double & sigma,
														unsigned int c)
{
	unsigned int n_pixel, n_pixel_mask;
	if ( width_step == 0)
		width_step = width;
	if ( mask_width_step == 0)
		mask_width_step = width;
		
	memset ( _data, 0, sizeof( double ) * _nb_bins ); // 0
	_nb_pixels = 0; // Nb pixels
	
	//Calcul des valeurs qui bornent l'histogramme
	for ( unsigned int i = 0; i < height; ++ i )
	{
		for ( unsigned int j = 0; j < width; ++ j )
		{
			n_pixel = i * width_step + j;
			n_pixel_mask = i * mask_width_step + j;
			if (  data[ n_pixel ] < _nb_bins && mask_data[ n_pixel_mask] != 0 )
			{

				_nb_pixels ++;
				_data[ (unsigned char) data[ n_pixel ] ]++;
			}
		}
	}
	if ( sigma != 0 )
	{
		filter (	c,
					sigma );	
		compute_derivates ();
	}
}


template<class type1, class type2> int c_histogram :: correct_alpha(	type2 * data_out,
																		const type1 * data_in,
																		unsigned int width,
																		unsigned int height,
																		unsigned int width_step_in,
																		unsigned int width_step_out,
																		double new_mean )
{
	//Calcul de la moyenne
	double m = 0;
	unsigned int nb_pixels = 0;
	for ( unsigned int i = 0; i < _nb_bins; ++ i )
	{
		m += _data[i] * i;
		nb_pixels += _data[i];
	}
	
	m /= nb_pixels;
	
	new_mean *= (_nb_bins - 1);
	double max_mean = (_nb_bins - 1) * ( nb_pixels - _data[0]) / ( (double) nb_pixels ),
		   min_mean = (_nb_bins - 1) * _data[_nb_bins - 1] / ( (double) nb_pixels );
	
	

	if ( max_mean <= new_mean )
	{
		for ( unsigned int i = 0; i < height; ++ i )
		{
			for ( unsigned int k = 0; k < width; ++ k ) 
			{
				unsigned int n_pixel_in,
							 n_pixel_out;
				n_pixel_out = i * width_step_out + k;
				n_pixel_in = i * width_step_in + k;
				double v = data_in[ n_pixel_in ];
				if ( v == 0 )
					data_out[ n_pixel_out ] = (type1) 0;
				else
					data_out[ n_pixel_out ] = (type1) _nb_bins - 1;
				
			}
		}
		return -1;
	}
	
	if ( min_mean >= new_mean )
	{
		for ( unsigned int i = 0; i < height; ++ i )
		{
			for ( unsigned int k = 0; k < width; ++ k ) 
			{
				unsigned int n_pixel_in,
							 n_pixel_out;
				n_pixel_out = i * width_step_out + k;
				n_pixel_in = i * width_step_in + k;
				double v = data_in[ n_pixel_in ];
				if ( v == (_nb_bins - 1) )
					data_out[ n_pixel_out ] = (type1) _nb_bins - 1;
				else
					data_out[ n_pixel_out ] = (type1) 0;
			}
		}
		return 1;
	}
	


	
	type2 * new_values = new type2[_nb_bins];
	double alpha;
	double mean = m,
		   previous_mean = 0;
	
	double alpha_min, 
		   alpha_max;
	if ( new_mean > m )
	{
		alpha_min = 0;
		alpha_max = 1;
		
		do
		{
			alpha = ( alpha_max + alpha_min ) / 2;
			//Calcul des nouvelles valeurs
			for ( unsigned int i = 0; i < _nb_bins; ++ i )
				new_values[i] = (type2) (_nb_bins - 1 ) * pow( i / ( (double) (_nb_bins - 1) ), alpha );
			
			//Calcul de la nouvelle moyenne
			previous_mean = mean;
			
			mean = 0;
			for ( unsigned int i = 0; i < _nb_bins; ++ i )
				mean += new_values[i] * _data[i];
			mean /= _nb_pixels;
			
			if ( mean > new_mean )
				alpha_min = alpha;
			else
				alpha_max = alpha;
			
		} while ( abs( mean - new_mean ) > 1 && abs(mean - previous_mean) > 1);
	}
	else 
	{
		alpha_min = 0;
		alpha_max = 1;
		
		do
		{
			alpha = ( alpha_max + alpha_min ) / 2;
			//Calcul des nouvelles valeurs
			for ( unsigned int i = 0; i < _nb_bins; ++ i )
				new_values[i] = (type2) ( _nb_bins - 1 )* pow( i / ( (double) (_nb_bins - 1 ) ), 1.0 / alpha );
			
			//Calcul de la nouvelle moyenne
			previous_mean = mean;
			
			mean = 0;
			for ( unsigned int i = 0; i < _nb_bins; ++ i )
				mean += new_values[i] * _data[i];
			mean /= _nb_pixels;
			
			if ( mean > new_mean )
				alpha_max = alpha;
			else
				alpha_min = alpha;
			
		} while ( abs( mean - new_mean ) > 1 && abs(mean - previous_mean) > 1);
		alpha = 1/alpha;
	}
	
	for ( unsigned int i = 0; i < height; ++ i )
	{
		for ( unsigned int k = 0; k < width; ++ k ) 
		{
			unsigned int n_pixel_in,
						 n_pixel_out;
			n_pixel_out = i * width_step_out + k;
			n_pixel_in = i * width_step_in + k;
			double v = data_in[ n_pixel_in ];
			data_out[ n_pixel_out ] = (type1) new_values[ (unsigned char) v];
			
		}
	}
	
	delete[] new_values;
	return 0;
}

template<class type> int c_histogram :: correct_alpha(	type * data,
														unsigned int width,
														unsigned int height,
														unsigned int width_step,
														double new_mean )
{
	return correct_alpha ( 	data, 
							(const type*) data, 
							width, 
							height, 
							width_step, 
							width_step, 
							new_mean );
	
	
	
	
}









template <class type1, class type2> void c_histogram :: equalize(	type2 * data_out,
															const type1 * data_in,
															unsigned int width,
															unsigned int height,
															unsigned int width_step_in,
															unsigned int width_step_out)
{
	unsigned char * thresholds = new unsigned char[_nb_bins],
				  * toto = new unsigned char[_nb_bins];
	memset( thresholds, 0, sizeof(char) * _nb_bins );
	memset( toto, 0, sizeof(char) * _nb_bins );
	
	double 	s = ( (double) _nb_pixels ) / _nb_bins,
			cum = 0;
	unsigned int j =  0;
	for ( unsigned int i = 0; i < _nb_bins; ++ i )
	{
		cum += _data[i];
		while ( cum >= (j + 1) * s )
		{
			thresholds[j] = i;
			++ j;
			
		}


	}
	thresholds[_nb_bins - 1 ] = 255;
	j = 0;
	for ( unsigned int i = 0; i < _nb_bins; ++ i )
	{
		for (; j <= thresholds[i]; ++ j)
		{
			toto[j] = i;
		}

	}
	
	for ( unsigned int i = 0; i < height; ++ i )
	{
		for ( unsigned int k = 0; k < width; ++ k ) 
		{
			unsigned int n_pixel_in,
						 n_pixel_out;
			n_pixel_out = i * width_step_out + k;
			n_pixel_in = i * width_step_in + k;
			double v = data_in[ n_pixel_in ];
			data_out[ n_pixel_out ] = (type1) toto[ (unsigned char) v];
			
		}
	}
	
	delete[] thresholds;
	delete[] toto;
}














template <class type> void c_histogram :: equalize(	type * data,
													unsigned int width,
													unsigned int height,
													unsigned int width_step)
{
	equalize( data, (const type*) data, width, height, width_step, width_step );
	
}


void c_histogram :: filter (	unsigned int c,
								const double & sigma )
{
	double * gaussian_kernel;
	double variance_inv;
	unsigned int m;
	double threshold = _nb_pixels / (10.0 * _nb_bins);
	
	if ( c == 0)
		c = 3;
	//Taille de la fenêtre gaussienne
	m = (unsigned int) c * (sigma) + 1; // M
	
	//Noyeau gaussian
	gaussian_kernel = new double[ 2 * m + 1 ]; //n=1 to M
	variance_inv = 1.0 / ( sigma * sigma);
	gaussian_kernel[m] = 1;
	for (unsigned int n = 1; n <= m; ++ n)
	{
		gaussian_kernel[m + n] = exp( -(double)  (n * n) * variance_inv );
		gaussian_kernel[m - n] = exp( -(double)  (n * n) * variance_inv );
	}

	filter_1d( _s_data, _data, _nb_bins, gaussian_kernel, 2 * m + 1 );
	//Suppression des zones avec une valeur trop faible
	for ( unsigned int i = 0; i < _nb_bins; ++ i)
	{
		if ( _s_data[i] < threshold )
		{
			_s_data[i] = 0;
		}
	}
	
	delete[] gaussian_kernel;
	
}


template <class type> int c_histogram :: stretch( 	type * data,
													unsigned int width,
													unsigned int height,
													unsigned int width_step,
													double eps_1,
													double eps_2 )
{
	return stretch( 	data, 
						(const type*) data, 
						width, 
						height, 
						width_step, 
						width_step, 
						eps_1,
						eps_2 );
}

template<class type1, class type2> int c_histogram :: stretch( 	type1 * data,
																const type2 * data_in,
																unsigned int width,
																unsigned int height,
																unsigned int width_step_in,
																unsigned int width_step_out,
																double eps_low,
																double eps_up )
{
	unsigned int nb_pixel = 0;
	unsigned int min = 0, max;
	
	double n_lost_pixels_low,
		   n_lost_pixels_up;
	
	
	if ( width_step_in == 0 )
		width_step_in = width;
	if ( width_step_out == 0 )
		width_step_out = width;
		
	if ( eps_low < 0 || eps_low >= 1.0) //Perte de 5% max
		return 1;
	if ( eps_up < 0 || eps_up >= 1.0) //Perte de 5% max
		return 1;
		
	n_lost_pixels_low = eps_low * _nb_pixels;
	n_lost_pixels_up = eps_up * _nb_pixels;
	
	for ( min = 0; min < _nb_bins; ++ min )
	{
		nb_pixel += _data[min];
		if ( nb_pixel > n_lost_pixels_low )
			break;
	}
		
	nb_pixel = 0;
	for ( max = _nb_bins - 1; max != min; -- max )
	{
		nb_pixel += _data[max];
		if ( nb_pixel > n_lost_pixels_up )
			break;
	}
	if ( max - min < 32 ) 
	{
		memset ( 	data, 
					0, 
					width_step_out * sizeof( type1) * height );
		return 0;
		
	}
	//Coeff
	double a, b;
	a = 255.0 / ( (double) (max - min) );
	b = - (min * a);

	//Mise à niveau de l'histogramme
	for ( unsigned int i = 0; i < height; ++ i )
	{
		for ( unsigned int j = 0; j < width; ++ j ) 
		{
			unsigned int n_pixel_in,
						 n_pixel_out;
			double v;
			n_pixel_out = i * width_step_out + j;
			n_pixel_in = i * width_step_in + j;
			v = a * data_in[ n_pixel_in ] + b;
			if ( v > 255 )
				v = 255;
			else if ( v < 0 )
				v = 0;
			data[ n_pixel_out ] = (type1) ( v );
		}
		
	}
	return 0;
}

void c_histogram :: compute_derivates ()
{
	unsigned int end = _nb_bins - 1;
	double threshold = _nb_pixels / (100.0 * _nb_bins);
	//D1
	for ( unsigned int i = 0; i < end; ++ i)
	{
		_d_data[i] = ( _s_data[i + 1] - _s_data[i] );
		if (abs(_d_data[i]) < threshold)
			_d_data[i] = 0;
	}

	//D2
	end --;
	//0
	for ( unsigned int i = 0; i < end; ++ i)
	{
		_d_d_data[i] = ( _d_data[i + 1] - _d_data[i] );
		if (abs(_d_d_data[i]) < threshold)
			_d_d_data[i] = 0;
		
	}
	


	
}

void c_histogram :: search_maxima (	unsigned int * maxima,
									unsigned int & nb_min,
									unsigned int nb_min_max ) const
{
	nb_min = 0;
	if (nb_min_max == 0)
		return;

	int last_sign = _d_data[0];
	unsigned int end = _nb_bins - 1;
	unsigned int i = 1;
	
	
	//Gauche

	for (; i < end; ++ i )
	{
		if ( last_sign >= 0 )
		{
			if ( _d_data[i] < 0)
			{
				unsigned int j;
				if ( nb_min != nb_min_max )
				{
					for ( j = 0; j < nb_min; ++ j )
					{
						if ( _s_data[i] > _s_data[maxima[j]] )
						{
							break;
						}
					}
					for (unsigned int k = nb_min; k != j; -- k)
					{
						maxima[k] = maxima[k - 1];
					}
					maxima[j] = i;
					nb_min ++ ;
				}
				else
				{
					for ( j = 0; j < nb_min; ++ j )
					{
						if ( _s_data[i] > _s_data[maxima[j]] )
						{
							for (unsigned int k = nb_min - 1; k != j; -- k)
							{
								maxima[k] = maxima[k - 1];
							}
							maxima[j] = i;
							break;
						}
					}
				}
			}
			
		}
		last_sign = _d_data[i];
	}
	
	sort(maxima, maxima + nb_min);
}

void c_histogram :: search_minima (	unsigned int * minima,
									unsigned int & nb_min,
									unsigned int nb_min_max ) const
{
	nb_min = 0;
	if (nb_min_max == 0)
		return;
	nb_min = 1;
	minima[0] = 0;
	
	int last_sign = _d_data[1];
	unsigned int end = _nb_bins - 1;
	unsigned int i = 2;
	
	
	//Gauche

	for (; i < end - 1; ++ i )
	{
		int q = 0;
		if ( last_sign < 0 )
		{
			if ( _d_data[i] >= 0    )
			{
				q = 1;
			}
		}
		//~ if ( !q )
		//~ {
			//~ if ( _d_d_data[i + 1] * _d_d_data[i] < 0 )
			//~ {
				//~ q = 1;
			//~ }
		//~ } 
		if ( q )
		{
			unsigned int j;
			if ( nb_min != nb_min_max )
			{
				for ( j = 0; j < nb_min; ++ j )
				{
					if ( _s_data[i] < _s_data[minima[j]] )
					{
						break;
					}
				}
				for (unsigned int k = nb_min; k != j; -- k)
				{
					minima[k] = minima[k - 1];
				}
				minima[j] = i;
				nb_min ++ ;
			}
			else
			{
				for ( j = 0; j < nb_min; ++ j )
				{
					if ( _s_data[i] < _s_data[minima[j]] )
					{
						for (unsigned int k = nb_min - 1; k != j; -- k)
						{
							minima[k] = minima[k - 1];
						}
						minima[j] = i;
						break;
					}
				}
			}
		
		}
		last_sign = _d_data[i];
	}
	
	sort(minima, minima + nb_min);
}

c_histogram :: ~c_histogram()
{
	free();
	initialize();
}

void c_histogram :: initialize()
{
	_nb_pixels = 0;
	_nb_bins = 0;
	_data = 0;
	_d_d_data = 0;
	_d_data = 0;
	_s_data = 0;
}

void c_histogram :: free()
{
	if ( _data )
		delete[] _data;
	if ( _s_data )
		delete[] _s_data;
	if ( _d_data )
		delete[] _d_data;
	if ( _d_d_data)
		delete[] _d_d_data;
}

template <class type> void c_histogram :: compute(	const type * data,
													unsigned int width,
													unsigned int height,
													unsigned int width_step,
													const double & sigma,
													unsigned int c)
{
	compute_histogram (	data,
						width,
						height,
						width_step);
	if ( sigma != 0 )
	{
		filter (	c,
					sigma );	
		compute_derivates ();
	}
}

/**@fn
 * @brief
 * Retourne la moyenne
 * 
 **/
double c_histogram :: get_mean( ) const
{
	double m = 0;
	if ( _nb_pixels == 0 )
		return -1;
		
	for ( unsigned int i = 0; i <_nb_bins; ++ i )
	{
		m += _data[i] * i;
	}
	
	m /= _nb_pixels;
	return m;
	
}


/**@fn
 * @brief
 * Retourne la moyenne
 * 
 **/
double c_histogram :: get_variance( ) const
{
	double m = get_mean( );
	if ( m == -1)
		return -1;
	double v = 0;
		
	for ( unsigned int i = 0; i <_nb_bins; ++ i )
	{
		v += _data[i] * pow( i - m, 2 );
	}
	
	m /= _nb_pixels;
	return m;
}



double c_histogram :: get_percentile( unsigned int n ) const
{
	double r = (n / 100.0) * _nb_pixels;
	if ( _nb_pixels == 0 || n > 100)
		return -1;
	double cum = 0;
	unsigned int i = 0;
	for ( ; i < _nb_bins; ++ i )
	{
		cum += _data[i];
		if ( cum >= r )
			break;
	}
	return i;  
}


C_HISTOGRAM_STRETCH_1( unsigned char );
C_HISTOGRAM_STRETCH_1( float );
C_HISTOGRAM_STRETCH_1( double );  

C_HISTOGRAM_COMPUTE( unsigned char )
C_HISTOGRAM_COMPUTE( float )
C_HISTOGRAM_COMPUTE( double )

C_HISTOGRAM_COMPUTE_HIST( unsigned char )
C_HISTOGRAM_COMPUTE_HIST( float )
C_HISTOGRAM_COMPUTE_HIST( double )

C_HISTOGRAM_STRETCH_2( unsigned char, unsigned char );
C_HISTOGRAM_STRETCH_2( unsigned char, float );
C_HISTOGRAM_STRETCH_2( unsigned char, double );
C_HISTOGRAM_STRETCH_2( float, unsigned char );
C_HISTOGRAM_STRETCH_2( float, float );
C_HISTOGRAM_STRETCH_2( float, double );
C_HISTOGRAM_STRETCH_2( double, unsigned char );
C_HISTOGRAM_STRETCH_2( double, float );
C_HISTOGRAM_STRETCH_2( double, double );

C_HISTOGRAM_EQ( unsigned char )
C_HISTOGRAM_EQ( float )
C_HISTOGRAM_EQ( double )

C_HISTOGRAM_EQ_2( unsigned char, unsigned char )
C_HISTOGRAM_EQ_2( unsigned char, float )
C_HISTOGRAM_EQ_2( unsigned char, double )

C_HISTOGRAM_EQ_2( float, unsigned char )
C_HISTOGRAM_EQ_2( float, float )
C_HISTOGRAM_EQ_2( float, double )

C_HISTOGRAM_EQ_2( double, unsigned char )
C_HISTOGRAM_EQ_2( double, float )
C_HISTOGRAM_EQ_2( double, double )

C_HISTOGRAM_C_ALPHA( unsigned char )
C_HISTOGRAM_C_ALPHA( float )
C_HISTOGRAM_C_ALPHA( double )

C_HISTOGRAM_C_ALPHA_2( unsigned char, unsigned char )
C_HISTOGRAM_C_ALPHA_2( unsigned char, float )
C_HISTOGRAM_C_ALPHA_2( unsigned char, double )

C_HISTOGRAM_C_ALPHA_2( float, unsigned char )
C_HISTOGRAM_C_ALPHA_2( float, float )
C_HISTOGRAM_C_ALPHA_2( float, double )

C_HISTOGRAM_C_ALPHA_2( double, unsigned char )
C_HISTOGRAM_C_ALPHA_2( double, float )
C_HISTOGRAM_C_ALPHA_2( double, double )
