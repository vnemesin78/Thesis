#include "tkalman_nc_em_base.hpp"
tkalman_nc_em_base :: tkalman_nc_em_base(	const gsl_vector * t0,
											const gsl_matrix * sqrt_q0,
											const gsl_matrix * f,
											const gsl_matrix * sqrt_q,
											unsigned int size_x,
											unsigned int n_max,
											unsigned int nb_signal_max,
											const gsl_vector * x_mask,
											const gsl_vector * y_mask,
											bool estimate_initial_state ) throw (exception &)
{
	initialize();
	//Validité des arguments d'entrée
	if ( !(t0 && sqrt_q0 && size_x && n_max && nb_signal_max) )
		throw (invalid_argument("t0, sqrt_q0, f, sqrt_q, size_x, p_max or n_max are NULL.\n"));
	//Dim
	if ( (sqrt_q0->size1 != sqrt_q0->size2))
		throw (invalid_argument("sqrt_q0, f or sqrt_q are not square matrix.\n"));
		
	if ( (sqrt_q0->size1 != t0->size ) )
		throw (invalid_argument("dim(q0), dim(Q), dim(F) or dim(t0) are different.\n"));
	
	unsigned int size_t, 
				 size_y;
	size_t = sqrt_q0->size2;
	size_y = size_t - size_x;
	
	//Vérif de la dim. de y;
	if (!size_y or !size_x)
		throw(invalid_argument("size_y or size_x are 0!\n"));
	

	free();
	initialize();
	_size_x = size_x;
	_size_y = size_y;
	_size_t = size_t;
	_n_max = n_max;
	_nb_signal_max = nb_signal_max;
	_estimate_initial_state = estimate_initial_state;
	try
	{
		alloc();
	}
	catch (exception & except)
	{
		throw(except);
	}
	//Recopie des paramètres
	{
		gsl_vector_memcpy(_t0, t0);
		gsl_matrix_memcpy(_sqrt_q0, sqrt_q0);
	}
	//Création du masque
	{
		if (x_mask == NULL)
			gsl_vector_set_all(_x_mask, 1);
		else
		{
			if (x_mask->size != _size_t)
			{
				throw (invalid_argument("dim of x_mask is not _size_t.\n"));
			}
			gsl_vector_memcpy(_x_mask, x_mask);
		}
		if (y_mask == NULL)
			gsl_vector_set_all(_y_mask, 1);
		else
		{
			if (y_mask->size != _size_t)
			{
				throw (invalid_argument("dim of y_mask is not _size_t.\n"));
			}
			gsl_vector_memcpy(_y_mask, y_mask);
		}
	}
	//Création des parties constantes de F
	{
		{
			gsl_matrix_const_view v1 = gsl_matrix_const_submatrix(f, 0, 0, _size_x, _size_t);
			gsl_matrix_transpose_memcpy(_f_xt_t, &(v1.matrix));
		}
		{
			gsl_matrix_const_view v1 = gsl_matrix_const_submatrix(f, _size_x, 0, _size_y, _size_t);
			gsl_matrix_transpose_memcpy(_f_yt_t, &(v1.matrix));
		}
	}
	
	//Objets
	try
	{
		create_object();
	}
	catch (exception & except)
	{
		throw(except);
	}
}



void tkalman_nc_em_base :: setup (	const gsl_vector * t0,
									const gsl_matrix * sqrt_q0,
									const gsl_matrix * f,
									const gsl_matrix * sqrt_q,
									unsigned int size_x,
									unsigned int n_max,
									unsigned int nb_signal_max,
									const gsl_vector * x_mask,
									const gsl_vector * y_mask,
									bool estimate_initial_state ) throw (exception &)
{
	//Validité des arguments d'entrée
	if ( !(t0 && sqrt_q0 && size_x && n_max && nb_signal_max) )
		throw (invalid_argument("t0, sqrt_q0, f, sqrt_q, size_x, p_max or n_max are NULL.\n"));
	//Dim
	if ( (sqrt_q0->size1 != sqrt_q0->size2))
		throw (invalid_argument("sqrt_q0, f or sqrt_q are not square matrix.\n"));
		
	if ( (sqrt_q0->size1 != t0->size ) )
		throw (invalid_argument("dim(q0), dim(Q), dim(F) or dim(t0) are different.\n"));
	
	unsigned int size_t, 
				 size_y;
	size_t = sqrt_q0->size2;
	size_y = size_t - size_x;
	
	//Vérif de la dim. de y;
	if (!size_y or !size_x)
		throw(invalid_argument("size_y or size_x are 0!\n"));
	

	free();
	initialize();
	_size_x = size_x;
	_size_y = size_y;
	_size_t = size_t;
	_n_max = n_max;
	_nb_signal_max = nb_signal_max;
	_estimate_initial_state = estimate_initial_state;
	try
	{
		alloc();
	}
	catch (exception & except)
	{
		throw(except);
	}
	//Recopie des paramètres
	{
		gsl_vector_memcpy(_t0, t0);
		gsl_matrix_memcpy(_sqrt_q0, sqrt_q0);
	}
	//Création du masque
	{
		if (x_mask == NULL)
			gsl_vector_set_all(_x_mask, 1);
		else
		{
			if (x_mask->size != _size_t)
			{
				throw (invalid_argument("dim of x_mask is not _size_t.\n"));
			}
			gsl_vector_memcpy(_x_mask, x_mask);
		}
		if (y_mask == NULL)
			gsl_vector_set_all(_y_mask, 1);
		else
		{
			if (y_mask->size != _size_t)
			{
				throw (invalid_argument("dim of y_mask is not _size_t.\n"));
			}
			gsl_vector_memcpy(_y_mask, y_mask);
		}
	}
	//Création des parties constantes de F
	{
		{
			gsl_matrix_const_view v1 = gsl_matrix_const_submatrix(f, 0, 0, _size_x, _size_t);
			gsl_matrix_transpose_memcpy(_f_xt_t, &(v1.matrix));
		}
		{
			gsl_matrix_const_view v1 = gsl_matrix_const_submatrix(f, _size_x, 0, _size_y, _size_t);
			gsl_matrix_transpose_memcpy(_f_yt_t, &(v1.matrix));
		}
	}
	
	//Objets
	try
	{
		create_object();
	}
	catch (exception & except)
	{
		throw(except);
	}
}


tkalman_nc_em_base :: ~tkalman_nc_em_base()
{
	free();
	initialize();
}

void tkalman_nc_em_base :: initialize()
{
	f_tools_x = 0;
	f_tools_y = 0;
	//Fusion des résultats
	fusion = 0;
		
	//Sommes
	_nb_signal_max = 0;
		
	//Paramètres
	_size_x = 0;
	_size_y = 0;
	_size_t = 0;
	_n_max = 0;
	_t0 = 0;
	_sqrt_q0 = 0;
		
	_x_mask = 0;
	_y_mask = 0;
	_estimate_initial_state = false;
	//Tmp
	tmp_f_xt_t = 0;
	tmp_f_yt_t = 0;
	_f_xt_t = 0;
	_f_yt_t = 0;
	t_0s = 0;
	sqrt_q_0s = 0;
}

void tkalman_nc_em_base :: alloc() throw(exception &)
{
	if (!_t0)
	{
		try
		{
			_t0 = gsl_vector_alloc(_size_t);
		}
		catch(exception & e)
		{
			throw (e);
		}
	}
	if (!_f_xt_t)
	{
		try
		{
			_f_xt_t = gsl_matrix_alloc(_size_t, _size_x);
		}
		catch(exception & e)
		{
			throw (e);
		}
	}
	
	if (!_f_yt_t)
	{
		try
		{
			_f_yt_t = gsl_matrix_alloc(_size_t, _size_y);
		}
		catch(exception & e)
		{
			throw (e);
		}
	}
	if (!_x_mask)
	{
		try
		{
			_x_mask = gsl_vector_alloc(_size_t);
		}
		catch(exception & e)
		{
			throw (e);
		}
	}
	
	if (!_y_mask)
	{
		try
		{
			_y_mask = gsl_vector_alloc(_size_t);
		}
		catch(exception & e)
		{
			throw (e);
		}
	}
	
	
	if (!_sqrt_q0)
	{
		try
		{
			_sqrt_q0 = gsl_matrix_alloc(_size_t, _size_t);
		}
		catch(exception & e)
		{
			throw (e);
		}
	}
	
	if (!tmp_f_xt_t)
	{
		try
		{
			tmp_f_xt_t = gsl_matrix_alloc(_size_t, _size_x);
		}
		catch(exception & e)
		{
			throw (e);
		}
	}
	
	if (!tmp_f_yt_t)
	{
		try
		{
			tmp_f_yt_t = gsl_matrix_alloc(_size_t, _size_y);
		}
		catch(exception & e)
		{
			throw (e);
		}
	}

	if (!t_0s)
	{
		t_0s = new gsl_vector *[_nb_signal_max];
		for ( unsigned int i = 0; i < _nb_signal_max; ++ i )
		{
			t_0s[i] = gsl_vector_alloc( _size_t );
		}
	}
	
	if (!sqrt_q_0s)
	{
		sqrt_q_0s = new gsl_matrix *[_nb_signal_max];
		for ( unsigned int i = 0; i < _nb_signal_max; ++ i )
		{
			sqrt_q_0s[i] = gsl_matrix_alloc( _size_t, _size_t );
		}
	}
	

}

void tkalman_nc_em_base :: free()
{
	if ( t_0s )
	{
		for ( unsigned int i = 0; i < _nb_signal_max; ++ i )
		{
			gsl_vector_free( t_0s[i] );
		}
		delete[] t_0s;
	}
	
	if ( sqrt_q_0s )
	{
		for ( unsigned int i = 0; i < _nb_signal_max; ++ i )
		{
			gsl_matrix_free( sqrt_q_0s[i] );
		}
		delete[] sqrt_q_0s;
	}
	
	
	
	if (f_tools_x)
		delete f_tools_x;
	if (f_tools_y)
		delete f_tools_y;
	if (_f_xt_t)
	{
		gsl_matrix_free(_f_xt_t);
	}
	
	if (_f_yt_t)
	{
		gsl_matrix_free(_f_yt_t);
	}
	if (_t0)
	{
		gsl_vector_free(_t0);
	}
	
	if (_x_mask)
	{
		gsl_vector_free(_x_mask);
	}
	
	if (_y_mask)
	{
		gsl_vector_free(_y_mask);
	}
	
	if (_sqrt_q0)
	{
		gsl_matrix_free(_sqrt_q0);
	}
	
	if (tmp_f_xt_t)
	{
		gsl_matrix_free(tmp_f_xt_t);
	}
	
	if (tmp_f_yt_t)
	{
		gsl_matrix_free(tmp_f_yt_t);
	}
	
	
}


void tkalman_nc_em_base :: create_object() throw(exception &)
{
	try
	{
		f_tools_x = new auxi_function_tools(_x_mask, _f_xt_t);
	}
	catch (exception & e)
	{
		throw(e);
	}
	try
	{
		f_tools_y = new auxi_function_tools(_y_mask, _f_yt_t);
	}
	catch (exception & e)
	{
		throw(e);
	}
	try
	{
		fusion = new tkalman_nc_fusion( _size_x, _size_y );
	}
	catch(exception & e)
	{
		throw (e);
	}
}



