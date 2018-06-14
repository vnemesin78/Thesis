#include "tkalman_nc_filter_base.hpp"

tkalman_nc_filter_base :: tkalman_nc_filter_base ( 	const gsl_vector * t0,
													const gsl_matrix * sqrt_q0,
													unsigned int size_x,
													unsigned int n_max )  throw(exception &)
{
	initialize();
	//Validité des arguments d'entrée
	if ( !(t0 && sqrt_q0  && size_x && n_max) )
		throw (invalid_argument("t0, sqrt_q0, size_x or n_max are NULL.\n"));
	//Dim
	if ( (sqrt_q0->size1 != sqrt_q0->size2))
		throw (invalid_argument("sqrt_q0 is not square matrix.\n"));
		
	if ( sqrt_q0->size1 != t0->size  )
		throw (invalid_argument("dim(q0) and dim(t0) are different.\n"));
	unsigned int size_t, 
				 size_y;
	size_t = sqrt_q0->size2;
	size_y = size_t - size_x;
	
	//Vérif de la dim. de y;
	if (!size_y or !size_x)
		throw(invalid_argument("size_y or size_x are 0!\n"));
	

	free();
	initialize();
	_t0 = t0;
	_sqrt_q0 = sqrt_q0;
	
	_size_x = size_x;
	_size_y = size_y;
	_size_t = size_t;
	_n_max = n_max;
	
	alloc();

}

void tkalman_nc_filter_base :: set_n_max( unsigned int n )  throw(exception &)
{
	if (n != _n_max)
	{
		tkalman_covariance_unref(_sqrt_p_p, _n_max + 1);
		tkalman_covariance_unref(_sqrt_p_s, _n_max + 1);
		tkalman_covariance_unref(_sqrt_p_f, _n_max);
		tkalman_covariance_unref(_sqrt_q_pp, _n_max + 1);
		tkalman_covariance_unref(_c_s, _n_max);
		tkalman_covariance_unref(_sqrt_s, _n_max);
		
		tkalman_expectation_unref(_x_p, _n_max + 1);
		tkalman_expectation_unref(_x_s, _n_max + 1);
		tkalman_expectation_unref(_x_f, _n_max);
		tkalman_expectation_unref(_t_pp, _n_max + 1);
		tkalman_expectation_unref(_innovation, _n_max);
		
		_n_max = n;
		
		tkalman_expectation_ref_v2(_x_p, _n_max + 1, _size_x, _size_t);
		tkalman_expectation_ref_v2(_x_f, _n_max, _size_x, _size_t);
		tkalman_expectation_ref_v2(_x_s, _n_max + 1, _size_x, _size_t);
		tkalman_expectation_ref_v2(_t_pp, _n_max + 1, _size_t, _size_t);
		tkalman_expectation_ref_v2(_innovation, _n_max, _size_y, _size_y);
		
		
		tkalman_covariance_ref_v2(_sqrt_p_p, _n_max + 1, _size_x, _size_x, _size_t, _size_t);
		tkalman_covariance_ref_v2(_sqrt_p_f,_n_max, _size_x, _size_x, _size_t, _size_t);
		tkalman_covariance_ref_v2(_sqrt_p_s, _n_max + 1, _size_x, _size_x, _size_t, _size_t);
		tkalman_covariance_ref_v2(_sqrt_q_pp, _n_max + 1, _size_t, _size_t, _size_t, _size_t);
		tkalman_covariance_ref_v2(_c_s, _n_max, _size_x, _size_x, _size_x, _size_t);
		tkalman_covariance_ref_v2(_sqrt_s, _n_max, _size_y, _size_y, _size_y, _size_y);
	}
}


void tkalman_nc_filter_base :: setup (	const gsl_vector * t0,
										const gsl_matrix * sqrt_q0,
										unsigned int size_x,
										unsigned int n_max ) throw(exception &)
{
	free();
	initialize();
	
	//Validité des arguments d'entrée
	if ( !(t0 && sqrt_q0  && size_x && n_max) )
		throw (invalid_argument("t0, sqrt_q0, size_x or n_max are NULL.\n"));
	//Dim
	if ( (sqrt_q0->size1 != sqrt_q0->size2))
		throw (invalid_argument("sqrt_q0 is not square matrix.\n"));
		
	if ( sqrt_q0->size1 != t0->size  )
		throw (invalid_argument("dim(q0) and dim(t0) are different.\n"));
	
	unsigned int size_t, 
				 size_y;
	size_t = sqrt_q0->size2;
	size_y = size_t - size_x;
	
	//Vérif de la dim. de y;
	if (!size_y or !size_x)
		throw(invalid_argument("size_y or size_x are 0!\n"));
	
	_t0 = t0;
	_sqrt_q0 = sqrt_q0;
	
	_size_x = size_x;
	_size_y = size_y;
	_size_t = size_t;
	_n_max = n_max;
	
	alloc();
	
}

tkalman_nc_filter_base :: ~tkalman_nc_filter_base()
{
	free();
	initialize();
}


void tkalman_nc_filter_base :: free()
{
	tkalman_covariance_unref(_sqrt_p_p, _n_max + 1);
	tkalman_covariance_unref(_sqrt_p_s, _n_max + 1);
	tkalman_covariance_unref(_sqrt_p_f, _n_max);
	tkalman_covariance_unref(_sqrt_q_pp, _n_max + 1);
	tkalman_covariance_unref(_c_s, _n_max);
	tkalman_covariance_unref(_sqrt_s, _n_max);
	
	tkalman_expectation_unref(_x_p, _n_max + 1);
	tkalman_expectation_unref(_x_s, _n_max + 1);
	tkalman_expectation_unref(_x_f, _n_max);
	tkalman_expectation_unref(_t_pp, _n_max + 1);
	tkalman_expectation_unref(_innovation, _n_max);
		
	if (vect_t)
		gsl_vector_free(vect_t);
	if (mat_tt)
		gsl_matrix_free(mat_tt);
	if (prediction)
		delete prediction;
	if (filtering)
		delete filtering;
	if (smoothing)
		delete smoothing;
}

void tkalman_nc_filter_base :: initialize()
{
	_t0 = 0;
	_sqrt_q0 = 0;
    _size_x = 0;
    _size_y = 0;
    _size_t = 0;
    _n_max = 0;
    _t_pp = 0;
    _x_p = 0;
    _x_f = 0;
    _x_s = 0;
    _innovation = 0;
    _sqrt_q_pp = 0;
    _sqrt_p_p = 0;
	_sqrt_p_s = 0;
	_sqrt_p_f = 0;
	_sqrt_s = 0;
	_c_s = 0;
	prediction = 0;
	filtering = 0;
	smoothing = 0;
	vect_t = 0;
	mat_tt = 0;	
}

void tkalman_nc_filter_base ::  alloc() throw(exception &)
{
	try
	{
		unsigned int n = _n_max;
		_n_max = 0;
		set_n_max(n);
	}
	catch(exception & e)
	{
		throw(e);
	}
	
	if (!vect_t)
	{
		try
		{
			vect_t = gsl_vector_alloc(_size_t);
		}
		catch(exception & e)
		{
			throw(e);
		}
	}
	
	if (!mat_tt)
	{
		try
		{
			mat_tt = gsl_matrix_alloc(_size_t, _size_t);
		}
		catch(exception & e)
		{
			throw(e);
		}
	}

}


