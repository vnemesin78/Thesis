#include "fitting_ellipse.hpp"
#include <iostream> 
#include <cmath>
#include <cstring>
using namespace std;
fitting_ellipse :: fitting_ellipse ()
{
	initialize();
	m_matrix = gsl_matrix_alloc(6, 6);
	c_matrix = gsl_matrix_alloc(6, 6);
	e_vectors = gsl_matrix_complex_alloc(6,6);
	alpha = gsl_vector_complex_alloc(6);
	beta = gsl_vector_alloc(6);
	vector = gsl_vector_alloc(6);
	workspace = gsl_eigen_genv_alloc (6);
}

fitting_ellipse :: ~fitting_ellipse ()
{
	free();
	initialize();
}

void fitting_ellipse :: free()
{
	if ( m_matrix )
		gsl_matrix_free(m_matrix);
	if ( c_matrix)
		gsl_matrix_free(c_matrix);
	if (workspace)
		gsl_eigen_genv_free (workspace);
	if (e_vectors)
		gsl_matrix_complex_free(e_vectors);
	if (beta)
		gsl_vector_free(beta);
	if (alpha)
		gsl_vector_complex_free(alpha);
	if (vector)
		gsl_vector_free(vector);
}

void fitting_ellipse :: initialize()
{
	m_matrix = 0;
	c_matrix = 0;
	workspace = 0;
	m_x = 0;
	m_y = 0;
	s_x = 0;
	s_y = 0;
	beta = 0;
	alpha = 0;
	vector = 0;
	e_vectors = 0;
	memset( moments_4, 0, 5 * sizeof(double) );
	memset( moments_3, 0, 4 * sizeof(double) );
	memset( moments_2, 0, 3 * sizeof(double) );
	memset( moments_1, 0, 2 * sizeof(double) );
	memset( par, 0, 6 * sizeof(double) );
	_x_center = 0;
	_y_center = 0;
	_a = 0;
	_b = 0;
	_theta = 0;
}

#define FITTING_ELLIPSE_FIT(type) template void fitting_ellipse :: fit( const type * data,\
																		unsigned int n,\
																		unsigned int nb_iter,\
																		const unsigned char * mask);
template<class type> void fitting_ellipse :: fit( const type * data,
												  unsigned int n,
												  unsigned int nb_iter,
												  const unsigned char * mask)
{
	normalize_data(data, n, mask);
	for (unsigned int i = 0 ; i < nb_iter; ++ i)
	{
		compute_moments(data, n, mask);
		compute_singular_values();
		unnormalize();
		compute_ellipse_params();
	}
}

#define FITTING_ELLIPSE_DAT(type) template void fitting_ellipse :: normalize_data ( const type * data,\
																					unsigned int n,\
																				    const unsigned char * mask);
template<class type> void fitting_ellipse :: normalize_data ( const type * data,
															  unsigned int n,
															  const unsigned char * mask)
{
	unsigned int p = 0;
	m_x = 0;
	m_y = 0;
	
	//Moyenne
	for (unsigned int i = 0; i < n; ++ i)
	{
		int q = 0;
		if ( ! mask )
			q = 1;
		else if ( mask[i] )
			q = 1;
		if ( q )
		{
			m_x += data[2 * i];
			m_y += data[2 * i + 1];
			++ p;
		}
	}
	m_x /= p;
	m_y /= p;
	//Recherche du min et du max
	{
		double x_min = data[0],
			   x_max = data[0],
			   y_min = data[1],
			   y_max = data[1];
		
		for (unsigned int i = 1; i < n; ++ i)
		{
			int q = 0;
			if ( ! mask )
				q = 1;
			else if ( mask[i] )
				q = 1;
			if ( q )
			{
				if ( data[2 * i] > x_max )
					x_max = data[2 * i];
				else if ( data[2 * i] < x_min )
					x_min = data[2 * i];
					
				if ( data[2 * i + 1] > y_max )
					y_max = data[2 * i + 1];
				else if ( data[2 * i + 1] < y_min )
					y_min = data[2 * i + 1];
			} 
		}
		s_x = (x_max - x_min);
		s_y = (y_max - y_min);
	}
	
}

#define FITTING_ELLIPSE_MOMENTS(type) template void fitting_ellipse :: compute_moments( const type * data,\
																						unsigned int n,\
																					    const unsigned char * mask);
template<class type> void fitting_ellipse :: compute_moments( const type * data,
															  unsigned int n,
																const unsigned char * mask)
{
	unsigned int p = 0;
	memset( moments_4, 0, 5 * sizeof(double) );
	memset( moments_3, 0, 4 * sizeof(double) );
	memset( moments_2, 0, 3 * sizeof(double) );
	memset( moments_1, 0, 2 * sizeof(double) );
	
	for (unsigned int i = 0; i < n; ++ i)
	{
		int q = 0;
		if ( ! mask )
			q = 1;
		else if ( mask[i] )
			q = 1;
		
		if ( q )
		{
			++ p;
			double x1,
				   x2,
				   x3,
				   x4,
				   y1,
				   y2,
				   y3,
				   y4;
			x1 = ( data[2 * i] - m_x ) / s_x;
			x2 = x1 * x1;
			x3 = x2 * x1;
			x4 = x3 * x1;
			y1 =  ( data[2 * i + 1] - m_y ) / s_y;
			y2 = y1 * y1;
			y3 = y2 * y1;
			y4 = y3 * y1;
			moments_4[0] += x4;
			moments_4[1] += x3 * y1;
			moments_4[2] += x2 * y2;
			moments_4[3] += x1 * y3;
			moments_4[4] += 	 y4;
			
			moments_3[0] += x3;
			moments_3[1] += x2 * y1;
			moments_3[2] += x1 * y2;
			moments_3[3] += 	 y3;
			
			moments_2[0] += x2;
			moments_2[1] += x1 * y1;
			moments_2[2] +=      y2;
			
			moments_1[0] += x1;
			moments_1[1] += y1;
		}
	}
	//R 0
	m_matrix->data[0] = moments_4[0];	
	m_matrix->data[1] = moments_4[1];	
	m_matrix->data[2] = moments_4[2];	
	m_matrix->data[3] = moments_3[0];	
	m_matrix->data[4] = moments_3[1];	
	m_matrix->data[5] = moments_2[0];
	//R 1
	m_matrix->data[1 * m_matrix->tda] = moments_4[1];	
	m_matrix->data[1 * m_matrix->tda + 1] = moments_4[2];	
	m_matrix->data[1 * m_matrix->tda + 2] = moments_4[3];	
	m_matrix->data[1 * m_matrix->tda + 3] = moments_3[1];	
	m_matrix->data[1 * m_matrix->tda + 4] = moments_3[2];	
	m_matrix->data[1 * m_matrix->tda + 5] = moments_2[1];
	
	//R 2
	m_matrix->data[2 * m_matrix->tda] = moments_4[2];	
	m_matrix->data[2 * m_matrix->tda + 1] = moments_4[3];	
	m_matrix->data[2 * m_matrix->tda + 2] = moments_4[4];	
	m_matrix->data[2 * m_matrix->tda + 3] = moments_3[2];	
	m_matrix->data[2 * m_matrix->tda + 4] = moments_3[3];	
	m_matrix->data[2 * m_matrix->tda + 5] = moments_2[2];
	
	//R 3
	m_matrix->data[3 * m_matrix->tda] = moments_3[0];	
	m_matrix->data[3 * m_matrix->tda + 1] = moments_3[1];	
	m_matrix->data[3 * m_matrix->tda + 2] = moments_3[2];	
	m_matrix->data[3 * m_matrix->tda + 3] = moments_2[0];	
	m_matrix->data[3 * m_matrix->tda + 4] = moments_2[1];	
	m_matrix->data[3 * m_matrix->tda + 5] = moments_1[0];
	
	//R 4
	m_matrix->data[4 * m_matrix->tda] = moments_3[1];	
	m_matrix->data[4 * m_matrix->tda + 1] = moments_3[2];	
	m_matrix->data[4 * m_matrix->tda + 2] = moments_3[3];	
	m_matrix->data[4 * m_matrix->tda + 3] = moments_2[1];	
	m_matrix->data[4 * m_matrix->tda + 4] = moments_2[2];	
	m_matrix->data[4 * m_matrix->tda + 5] = moments_1[1];
	
	//R 5
	m_matrix->data[5 * m_matrix->tda] = moments_2[0];	
	m_matrix->data[5 * m_matrix->tda + 1] = moments_2[1];	
	m_matrix->data[5 * m_matrix->tda + 2] = moments_2[2];	
	m_matrix->data[5 * m_matrix->tda + 3] = moments_1[0];	
	m_matrix->data[5 * m_matrix->tda + 4] = moments_1[1];	
	m_matrix->data[5 * m_matrix->tda + 5] = p;
	
}

void fitting_ellipse :: compute_singular_values ()
{
	gsl_matrix_set_zero ( c_matrix );
	
	//Matrice de contrainte
	c_matrix->data[ 0 * c_matrix->tda + 2 ] = -2;
	c_matrix->data[ 1 * c_matrix->tda + 1 ] = 1;
	c_matrix->data[ 2 * c_matrix->tda + 0 ] = -2;
	
	//SVD
	
	gsl_eigen_genv ( m_matrix, 
					 c_matrix, 
					 alpha, 
					 beta, 
					 e_vectors, 
					 workspace );

	//Recup. des VP négatives
	for (unsigned int i = 0; i < 6; ++ i)
	{
		if ( beta->data[i] != 0 )
		{

			if ( alpha->data[i * alpha->stride] / beta->data[i] < 0 )
			{
				for (unsigned int j = 0; j < 6; ++ j)
				{
					vector->data[ j * vector->stride ] = e_vectors->data[2 * j * e_vectors->tda + i]; 

				}
				break;

			}
		}
	}

	
}


void fitting_ellipse :: unnormalize()
{
	par[0] = vector->data[0] * s_y * s_y;
	par[1] = vector->data[1 * vector->stride] * s_x * s_y;
	par[2] = vector->data[2 * vector->stride] * s_x * s_x;
	//-2*A(1)*sy*sy*mx - A(2)*sx*sy*my + A(4)*sx*sy*sy
	par[3] = -2 * vector->data[0] * s_y * s_y * m_x -  vector->data[1 * vector->stride] * s_x * s_y * m_y + vector->data[3 * vector->stride] * s_x * s_y * s_y;
	//-A(2)*sx*sy*mx - 2*A(3)*sx*sx*my + A(5)*sx*sx*sy
	par[4] = - vector->data[1 * vector->stride] * s_x * s_y * m_x - 2 * vector->data[2 * vector->stride] * s_x * s_x * m_y + vector->data[4 * vector->stride] * s_x * s_x * s_y;
	//A(1)*sy*sy*mx*mx + A(2)*sx*sy*mx*my + A(3)*sx*sx*my*my  - A(4)*sx*sy*sy*mx - A(5)*sx*sx*sy*my  + A(6)*sx*sx*sy*sy 
	par[5] = par[0] * m_x * m_x + par[1] * m_x * m_y + par[2] * m_y * m_y - m_y - vector->data[3 * vector->stride] * s_x * s_y * s_y * m_x - vector->data[4 * vector->stride] * s_x * s_x * s_y * m_y - vector->data[5 * vector->stride] * s_x * s_x * s_y * s_y;
	
}


void fitting_ellipse :: compute_ellipse_params()
{
	double cos_t,
		   sin_t,
		   sin2_t,
		   cos2_t,
		   sin_t_cos_t;
	
	double a_0,
		   a_x,
		   a_y,
		   a_xx,
		   a_yy;
	
	double tmp_x_center,
		   tmp_y_center,
		   tmp_r;
	_theta = 0.5 * atan2( par[1], par[0] - par[2] );
	cos_t = cos(_theta);
	sin_t = sin(_theta);
	sin2_t = sin_t * sin_t;
	cos2_t = cos_t * cos_t;
	sin_t_cos_t = sin_t * cos_t;
	

		   
	a_0 = par[5];
	a_x = par[3] * cos_t + par[4] * sin_t;
	a_y = - par[3] * sin_t + par[4] * cos_t;
	a_xx = par[0] * cos2_t + par[2] * sin2_t + par[1] * sin_t_cos_t;
	a_yy = par[0] * sin2_t + par[2] * cos2_t - par[1] * sin_t_cos_t;
	
	tmp_x_center = - a_x / ( 2 * a_xx );
	tmp_y_center = - a_y / ( 2 * a_yy );
	tmp_r = a_0 - a_xx * tmp_x_center * tmp_x_center - a_yy * tmp_y_center * tmp_y_center;
	
	_x_center = tmp_x_center * cos_t - tmp_y_center * sin_t;
	_y_center = tmp_x_center * sin_t + tmp_y_center * cos_t;
	m_x = _x_center;
	m_y = _y_center;
	_a = - tmp_r / (a_xx);
	_b = - tmp_r / (a_yy);
	if (_a > 0)
		_a = sqrt(_a);
	else
	{
		_a = sqrt(-_a);
	}
	if (_b > 0)
		_b = sqrt(_b);
	else
	{
		_b = sqrt(-_b);

	}
}

#define FITTING_ELLIPSE_CLASS(type) FITTING_ELLIPSE_FIT(type)\
									FITTING_ELLIPSE_DAT(type)\
									FITTING_ELLIPSE_MOMENTS(type)
									
FITTING_ELLIPSE_CLASS(double)
FITTING_ELLIPSE_CLASS(unsigned int)
